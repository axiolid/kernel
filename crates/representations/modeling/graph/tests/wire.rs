//! The geometry graph wire format (ADR 0085): round trips over every node
//! kind and every variant reachable from one, the golden payloads that
//! pin each published format version, and every refusal by name.

use std::collections::{BTreeMap, BTreeSet};
use std::path::PathBuf;

use axiolid_core::{
    Aabb, BooleanOperator, Frame2, Frame3, Interval, Plane3, Point3, Transform2, Transform3, Vec2,
    Vec3,
};
use axiolid_curve::{
    AngleGraph2, Axis, BSplineCurve, BSplineSurface, BankConvention, Banked3, Basis, Branch,
    CantForm, CantLaw, CantPiece, Carrier, Chain2, ChainPiece2, Circle2, Circle3, CurvatureLaw,
    Curve2, Curve3, Elevated3, ElevationLaw, Ellipse2, Ellipse3, Field2, Harmonic, ImplicitCell,
    ImplicitCurve2, ImplicitSection3, Intrinsic2, Intrinsic3, KnotSpec, LiftedCurve2, Line2, Line3,
    PairNode, PairSection3, PatchField2, Polyline, QuadraticGraph2, RailSide, RuledCarrier,
    RuledSection3, SeriesField2, Sinusoid2, TorusCarrier, TorusSection3, Trig2,
};
use axiolid_mesh::{AttributeChannel, Blend, NormalAttribute, PolygonFace, PolygonMesh, TriMesh};
use axiolid_model::wire::{self, FormatVersion, WireError, FORMAT_NAME, FORMAT_VERSION};
use axiolid_model::{
    BasisEnd, CurveRelation, CurveSegment, CurveStation, GeometryGraph, GeometryGraphBuilder,
    GeometryNode, GraphError, Instance, InstanceAtStation, MasterRepresentation, NodeId,
    OpenProfile, OrientedCurveStation, PointOnCurve, PointOnSurface, RunStation, SeamSide, Section,
    SectionAtStation, SolidOperation, Station, StationFrame, StationOffsets, StationOrientation,
    StationedOpenSection, StationedSection, SurfaceRelation, SurfaceSides, Transition,
    TrimSelector, TrimmingPreference,
};
use axiolid_primitive::{HalfSpace, Primitive};
use axiolid_profile::{
    CenterLineProfile, CircleProfile, Contour, ContourProfile, EllipseProfile, Profile,
    ProfileSegment, RectangleProfile, SectionProfile,
};
use axiolid_surface::{Cone, Cylinder, EllipticalCylinder, Plane, Sphere, Surface, Torus};
use axiolid_topology::{
    BRep, Edge, EdgeUse, Face, FaceBound, Loop, Orientation, Shell, Solid, Vertex,
};
use serde::de::DeserializeOwned;
use serde::Serialize;
use serde_json::{json, Value};

// ---------------------------------------------------------------------------
// Variant coverage
// ---------------------------------------------------------------------------

/// The wire tag of one enum value: its unit-variant string or the single
/// key of its map.
fn tag<T: Serialize>(value: &T) -> String {
    match serde_json::to_value(value).expect("a fixture value serialises") {
        Value::String(name) => name,
        Value::Object(map) if map.len() == 1 => map.keys().next().unwrap().clone(),
        other => panic!("not an externally tagged enum value: {other}"),
    }
}

/// Every variant name `T`'s deserialiser knows, read from the refusal of a
/// name it does not know. A variant added to `T` appears here whether or
/// not anyone remembered the wire.
fn variants_of<T: DeserializeOwned>() -> BTreeSet<String> {
    let message = serde_json::from_str::<T>("\"\\u0001probe\"")
        .err()
        .expect("the probe is no variant")
        .to_string();
    let expected = message
        .split_once(", expected ")
        .unwrap_or_else(|| panic!("not an unknown-variant message: {message}"))
        .1;
    expected
        .split('`')
        .skip(1)
        .step_by(2)
        .map(str::to_owned)
        .collect()
}

/// The variant tags a fixture actually wrote, per enum.
#[derive(Default)]
struct Coverage(BTreeMap<&'static str, BTreeSet<String>>);

impl Coverage {
    fn all<T: Serialize>(&mut self, name: &'static str, values: &[T]) {
        for value in values {
            self.one(name, value);
        }
    }

    fn one<T: Serialize>(&mut self, name: &'static str, value: &T) {
        self.0.entry(name).or_default().insert(tag(value));
    }
}

// ---------------------------------------------------------------------------
// The fixture: every node kind, every variant
// ---------------------------------------------------------------------------

fn frame2() -> Frame2 {
    Frame2 {
        origin: Vec2::new(0.5, -0.25),
        x: Vec2::X,
        y: Vec2::Y,
    }
}

fn frame3() -> Frame3 {
    Frame3 {
        origin: Vec3::new(1.0, 2.0, -3.0),
        x: Vec3::X,
        y: Vec3::Y,
        z: Vec3::Z,
    }
}

fn trig(seed: f64) -> Trig2 {
    Trig2 {
        constant: seed,
        cos: seed + 0.5,
        sin: -seed,
        cos2: 0.125,
        sin2: -0.0,
    }
}

fn knot_specs() -> Vec<KnotSpec> {
    vec![
        KnotSpec::Uniform,
        KnotSpec::QuasiUniform,
        KnotSpec::PiecewiseBezier,
        KnotSpec::Unspecified,
    ]
}

fn bspline2(knot_spec: KnotSpec) -> BSplineCurve<Vec2> {
    BSplineCurve {
        degree: 1,
        control_points: vec![Vec2::ZERO, Vec2::new(1.0, 0.5), Vec2::new(2.0, 0.0)],
        knots: vec![0.0, 1.0, 2.0],
        multiplicities: vec![2, 1, 2],
        weights: None,
        closed: false,
        self_intersect: Some(false),
        knot_spec,
    }
}

fn bspline_surface() -> BSplineSurface {
    BSplineSurface {
        u_degree: 1,
        v_degree: 1,
        control_points: vec![
            vec![Vec3::ZERO, Vec3::Y],
            vec![Vec3::X, Vec3::new(1.0, 1.0, 0.25)],
        ],
        u_knots: vec![0.0, 1.0],
        u_multiplicities: vec![2, 2],
        v_knots: vec![0.0, 1.0],
        v_multiplicities: vec![2, 2],
        weights: Some(vec![vec![1.0, 0.5], vec![0.5, 1.0]]),
        u_closed: false,
        v_closed: true,
        knot_spec: KnotSpec::Unspecified,
        self_intersect: None,
    }
}

fn curvature_laws() -> Vec<CurvatureLaw> {
    vec![
        CurvatureLaw::Constant { curvature: 0.1 },
        CurvatureLaw::Polynomial {
            coefficients: vec![0.0, 0.01],
        },
        CurvatureLaw::Sinusoid {
            mean: 0.0,
            amplitude: 0.02,
            angular_frequency: 0.5,
            phase: -0.0,
        },
        CurvatureLaw::Composite {
            polynomial: vec![0.001],
            harmonics: vec![Harmonic {
                amplitude: 0.01,
                angular_frequency: 2.0,
                phase: 0.25,
            }],
        },
        CurvatureLaw::Piecewise {
            breaks: vec![5.0],
            laws: vec![
                CurvatureLaw::Constant { curvature: 0.0 },
                CurvatureLaw::Polynomial {
                    coefficients: vec![0.0, 0.002],
                },
            ],
        },
    ]
}

fn elevation_laws() -> Vec<ElevationLaw> {
    vec![
        ElevationLaw::Polynomial {
            coefficients: vec![10.0, 0.02],
        },
        ElevationLaw::Piecewise {
            breaks: vec![50.0],
            laws: vec![
                ElevationLaw::Polynomial {
                    coefficients: vec![10.0, 0.02],
                },
                ElevationLaw::CircularArc {
                    height: 11.0,
                    grade: 0.02,
                    radius: -2000.0,
                },
            ],
        },
        ElevationLaw::CircularArc {
            height: 1.0,
            grade: -0.01,
            radius: 1500.0,
        },
        ElevationLaw::Intrinsic {
            height: 2.0,
            grade: 0.0,
            curvature: CurvatureLaw::Constant { curvature: 1e-4 },
        },
    ]
}

fn cant_forms() -> Vec<CantForm> {
    vec![
        CantForm::Polynomial {
            coefficients: vec![0.0, 0.1],
        },
        CantForm::Cosine {
            start: 0.1,
            change: 0.05,
        },
        CantForm::Sine {
            start: 0.15,
            change: -0.05,
        },
        CantForm::VienneseBend {
            start: 0.0,
            change: 0.07,
        },
        CantForm::AboutRail {
            rail: RailSide::Left,
            elevation: 0.0,
        },
        CantForm::AboutRail {
            rail: RailSide::Right,
            elevation: 0.5,
        },
    ]
}

fn ruled_carrier() -> RuledCarrier {
    RuledCarrier {
        frame: frame3(),
        x_radius: 2.0,
        y_radius: 1.5,
        slope: 0.0,
    }
}

fn torus_carrier() -> TorusCarrier {
    TorusCarrier {
        frame: frame3(),
        major_radius: 3.0,
        minor_radius: 1.0,
    }
}

fn carriers() -> Vec<Carrier> {
    vec![
        Carrier::Plane(frame3()),
        Carrier::Ruled(ruled_carrier()),
        Carrier::Sphere {
            frame: frame3(),
            radius: 2.0,
        },
        Carrier::Torus(torus_carrier()),
        Carrier::Spline(Box::new(bspline_surface())),
    ]
}

fn fields() -> Vec<Field2> {
    vec![
        Field2::Series(SeriesField2 {
            u: Basis::Power,
            v: Basis::Fourier,
            coefficients: vec![vec![1.0, 0.0, -1.0], vec![0.5]],
        }),
        Field2::Patches(PatchField2 {
            u_breaks: vec![0.0, 1.0],
            v_breaks: vec![0.0, 1.0],
            u_degree: 1,
            v_degree: 1,
            patches: vec![vec![1.0, -1.0, -1.0, 1.0]],
        }),
    ]
}

fn implicit_curve(field: Field2) -> ImplicitCurve2 {
    ImplicitCurve2 {
        field,
        cells: vec![
            ImplicitCell {
                axis: Axis::U,
                from: 0.0,
                to: 0.5,
                low: -1.0,
                high: 1.0,
                bridge: None,
            },
            ImplicitCell {
                axis: Axis::V,
                from: 0.5,
                to: 1.0,
                low: -1.0,
                high: 1.0,
                bridge: Some((0.25, 0.75)),
            },
        ],
    }
}

fn line2() -> Curve2 {
    Curve2::Line(Line2 {
        origin: Vec2::ZERO,
        direction: Vec2::X,
    })
}

fn line3() -> Curve3 {
    Curve3::Line(Line3 {
        origin: Vec3::ZERO,
        direction: Vec3::X,
    })
}

fn open_polyline2() -> Curve2 {
    Curve2::Polyline(Polyline {
        points: vec![Vec2::ZERO, Vec2::new(1.0, 0.5), Vec2::new(2.0, 0.0)],
        closed: false,
    })
}

fn closed_polyline2() -> Curve2 {
    Curve2::Polyline(Polyline {
        points: vec![Vec2::ZERO, Vec2::X, Vec2::new(1.0, 1.0)],
        closed: true,
    })
}

fn curve2s(coverage: &mut Coverage) -> Vec<Curve2> {
    let laws = curvature_laws();
    coverage.all("CurvatureLaw", &laws);
    let specs = knot_specs();
    coverage.all("KnotSpec", &specs);
    let fields = fields();
    coverage.all("Field2", &fields);
    coverage.all("Basis", &[Basis::Power, Basis::Fourier]);
    coverage.all("Axis", &[Axis::U, Axis::V]);
    let carriers = carriers();
    coverage.all("Carrier", &carriers);
    let pieces = vec![
        ChainPiece2::Intrinsic {
            curvature: CurvatureLaw::Constant { curvature: 0.0 },
            length: 10.0,
        },
        ChainPiece2::Parametric {
            curve: line2(),
            start: 0.0,
            length: 2.0,
        },
    ];
    coverage.all("ChainPiece2", &pieces);
    coverage.all("Branch", &[Branch::Plus, Branch::Minus]);

    let mut curves = vec![
        line2(),
        Curve2::Circle(Circle2 {
            frame: frame2(),
            radius: 1.25,
        }),
        Curve2::Ellipse(Ellipse2 {
            frame: frame2(),
            semi_axis_x: 2.0,
            semi_axis_y: 1.0,
        }),
        open_polyline2(),
        Curve2::Sinusoid(Sinusoid2 {
            mean: 1.0,
            cosine: 0.5,
            sine: -0.5,
        }),
        Curve2::QuadraticGraph(QuadraticGraph2 {
            a: trig(1.0),
            b: trig(0.0),
            c: trig(-1.0),
            branch: Branch::Plus,
        }),
        Curve2::AngleGraph(AngleGraph2 {
            a: trig(2.0),
            b: trig(1.0),
            c: trig(0.5),
            branch: Branch::Minus,
        }),
        Curve2::Chain(Chain2 {
            start: frame2(),
            pieces,
        }),
    ];
    curves.extend(
        specs
            .into_iter()
            .map(|spec| Curve2::BSpline(bspline2(spec))),
    );
    curves.extend(laws.into_iter().map(|curvature| {
        Curve2::Intrinsic(Intrinsic2 {
            start: frame2(),
            curvature,
            length: 20.0,
        })
    }));
    curves.extend(
        fields
            .into_iter()
            .map(|field| Curve2::Implicit(implicit_curve(field))),
    );
    curves.extend(carriers.into_iter().map(|carrier| {
        Curve2::Lifted(LiftedCurve2 {
            curve: Box::new(line3()),
            carrier,
            start: 0.0,
            end: 1.0,
            guide: vec![Vec2::ZERO, Vec2::new(0.5, 0.5)],
        })
    }));
    coverage.all("Curve2", &curves);
    curves
}

fn elevated(elevation: ElevationLaw) -> Elevated3 {
    Elevated3 {
        plan: Box::new(line2()),
        elevation,
    }
}

fn curve3s(coverage: &mut Coverage) -> Vec<Curve3> {
    let laws = elevation_laws();
    coverage.all("ElevationLaw", &laws);
    let forms = cant_forms();
    coverage.all("CantForm", &forms);
    coverage.all("RailSide", &[RailSide::Left, RailSide::Right]);
    let conventions = [
        BankConvention::TangentRotation,
        BankConvention::VerticalRise,
    ];
    coverage.all("BankConvention", &conventions);

    let mut curves = vec![
        line3(),
        Curve3::Circle(Circle3 {
            frame: frame3(),
            radius: 4.0,
        }),
        Curve3::Ellipse(Ellipse3 {
            frame: frame3(),
            semi_axis_x: 3.0,
            semi_axis_y: 2.0,
        }),
        Curve3::Polyline(Polyline {
            points: vec![Vec3::ZERO, Vec3::X, Vec3::new(1.0, 1.0, 1.0)],
            closed: false,
        }),
        Curve3::BSpline(BSplineCurve {
            degree: 1,
            control_points: vec![Vec3::ZERO, Vec3::Z],
            knots: vec![0.0, 1.0],
            multiplicities: vec![2, 2],
            weights: Some(vec![1.0, 2.0]),
            closed: false,
            self_intersect: None,
            knot_spec: KnotSpec::Unspecified,
        }),
        Curve3::Intrinsic(Intrinsic3 {
            start: frame3(),
            curvature: CurvatureLaw::Constant { curvature: 0.5 },
            torsion: CurvatureLaw::Constant { curvature: 0.1 },
            length: 6.0,
        }),
        Curve3::RuledSection(RuledSection3 {
            carrier: ruled_carrier(),
            graph: QuadraticGraph2 {
                a: trig(1.0),
                b: trig(2.0),
                c: trig(3.0),
                branch: Branch::Minus,
            },
        }),
        Curve3::TorusSection(TorusSection3 {
            torus: torus_carrier(),
            graph: AngleGraph2 {
                a: trig(0.25),
                b: trig(0.5),
                c: trig(0.75),
                branch: Branch::Plus,
            },
        }),
        Curve3::ImplicitSection(ImplicitSection3 {
            carrier: Carrier::Torus(torus_carrier()),
            curve: implicit_curve(fields().remove(0)),
        }),
        Curve3::PairSection(PairSection3 {
            first: Carrier::Spline(Box::new(bspline_surface())),
            second: Carrier::Plane(frame3()),
            nodes: vec![PairNode {
                point: Vec3::new(0.5, 0.5, 0.0),
                first: Vec2::new(0.5, 0.5),
                second: Vec2::new(0.5, 0.5),
            }],
        }),
    ];
    curves.extend(laws.into_iter().map(|law| Curve3::Elevated(elevated(law))));
    let (pivot, cant): (Vec<_>, Vec<_>) = forms
        .into_iter()
        .partition(|form| matches!(form, CantForm::AboutRail { .. }));
    let law = |forms: Vec<CantForm>| CantLaw {
        pieces: forms
            .into_iter()
            .map(|form| CantPiece { length: 25.0, form })
            .collect(),
    };
    let (cant, pivot) = (law(cant), law(pivot));
    curves.extend(conventions.into_iter().map(|convention| {
        Curve3::Banked(Banked3 {
            base: elevated(ElevationLaw::Polynomial {
                coefficients: vec![0.0],
            }),
            cant: cant.clone(),
            pivot: pivot.clone(),
            rail_head_distance: 1.5,
            convention,
        })
    }));
    coverage.all("Curve3", &curves);
    curves
}

fn surfaces(coverage: &mut Coverage) -> Vec<Surface> {
    let surfaces = vec![
        Surface::Plane(Plane { frame: frame3() }),
        Surface::Cylinder(Cylinder {
            frame: frame3(),
            radius: 1.0,
        }),
        Surface::EllipticalCylinder(EllipticalCylinder {
            frame: frame3(),
            semi_axis_x: 2.0,
            semi_axis_y: 1.0,
        }),
        Surface::Cone(Cone {
            frame: frame3(),
            radius: 1.0,
            semi_angle: 0.5,
        }),
        Surface::Sphere(Sphere {
            frame: frame3(),
            radius: 2.0,
        }),
        Surface::Torus(Torus {
            frame: frame3(),
            major_radius: 3.0,
            minor_radius: 0.5,
        }),
        Surface::BSpline(bspline_surface()),
    ];
    coverage.all("Surface", &surfaces);
    surfaces
}

fn section_profiles() -> Vec<SectionProfile> {
    vec![
        SectionProfile::I {
            depth: 0.3,
            width: 0.15,
            web_thickness: 0.007,
            flange_thickness: 0.011,
            fillet_radius: Some(0.015),
            flange_edge_radius: None,
            flange_slope: None,
        },
        SectionProfile::AsymmetricI {
            depth: 0.4,
            web_thickness: 0.01,
            bottom_flange_width: 0.2,
            bottom_flange_thickness: 0.015,
            bottom_fillet_radius: None,
            bottom_flange_edge_radius: None,
            bottom_flange_slope: None,
            top_flange_width: 0.15,
            top_flange_thickness: Some(0.012),
            top_fillet_radius: None,
            top_flange_edge_radius: None,
            top_flange_slope: Some(0.1),
        },
        SectionProfile::L {
            depth: 0.1,
            width: None,
            thickness: 0.01,
            fillet_radius: None,
            edge_radius: None,
            leg_slope: None,
        },
        SectionProfile::T {
            depth: 0.1,
            flange_width: 0.1,
            web_thickness: 0.01,
            flange_thickness: 0.01,
            fillet_radius: None,
            flange_edge_radius: None,
            web_edge_radius: None,
            web_slope: None,
            flange_slope: None,
        },
        SectionProfile::U {
            depth: 0.2,
            flange_width: 0.08,
            web_thickness: 0.006,
            flange_thickness: 0.01,
            fillet_radius: None,
            edge_radius: None,
            flange_slope: None,
        },
        SectionProfile::C {
            depth: 0.2,
            width: 0.07,
            wall_thickness: 0.003,
            girth: 0.02,
            internal_fillet_radius: None,
        },
        SectionProfile::Z {
            depth: 0.2,
            flange_width: 0.07,
            web_thickness: 0.005,
            flange_thickness: 0.008,
            fillet_radius: None,
            edge_radius: None,
        },
        SectionProfile::Trapezium {
            bottom_x: 1.0,
            top_x: 0.5,
            y: 0.75,
            top_offset: 0.25,
        },
    ]
}

fn contour() -> Contour {
    Contour {
        segments: vec![
            ProfileSegment {
                curve: line2(),
                domain: Interval {
                    start: 0.0,
                    end: 1.0,
                },
                same_sense: true,
            },
            ProfileSegment {
                curve: Curve2::Circle(Circle2 {
                    frame: frame2(),
                    radius: 0.5,
                }),
                domain: Interval {
                    start: 0.0,
                    end: core::f64::consts::PI,
                },
                same_sense: false,
            },
        ],
    }
}

fn rectangle() -> Profile {
    Profile::Rectangle(RectangleProfile {
        x: 2.0,
        y: 1.0,
        thickness: None,
        outer_radius: Some(0.1),
        inner_radius: None,
    })
}

fn profiles(coverage: &mut Coverage) -> Vec<Profile> {
    let sections = section_profiles();
    coverage.all("SectionProfile", &sections);
    let mut profiles = vec![
        rectangle(),
        Profile::Circle(CircleProfile {
            radius: 1.0,
            thickness: Some(0.1),
        }),
        Profile::Ellipse(EllipseProfile {
            semi_axis_x: 2.0,
            semi_axis_y: 1.0,
        }),
        Profile::Contour(ContourProfile {
            outer: contour(),
            holes: vec![contour()],
        }),
        Profile::Derived {
            basis: Box::new(rectangle()),
            transform: Transform2::from_cols_array(&[0.0, 1.0, -1.0, 0.0, 3.0, -4.0]),
        },
        Profile::Composite(vec![rectangle(), rectangle()]),
        Profile::CenterLine(CenterLineProfile {
            path: contour(),
            half_width: 0.1,
        }),
    ];
    profiles.extend(sections.into_iter().map(Profile::Section));
    coverage.all("Profile", &profiles);
    profiles
}

fn primitives(coverage: &mut Coverage) -> Vec<Primitive> {
    let primitives = vec![
        Primitive::Block {
            x: 1.0,
            y: 2.0,
            z: 3.0,
        },
        Primitive::Sphere { radius: 1.0 },
        Primitive::Cylinder {
            radius: 1.0,
            height: 2.0,
        },
        Primitive::Cone {
            radius: 1.0,
            height: 2.0,
        },
        Primitive::Pyramid {
            x: 1.0,
            y: 1.0,
            height: 2.0,
        },
        Primitive::Torus {
            major_radius: 3.0,
            minor_radius: 1.0,
        },
        Primitive::Wedge {
            x: 2.0,
            y: 1.0,
            height: 1.0,
            top_x_min: 0.0,
            top_x_max: 1.0,
            top_y_min: 0.0,
            top_y_max: 1.0,
        },
    ];
    coverage.all("Primitive", &primitives);
    primitives
}

/// A transform whose twelve entries are all distinct, so a column swap or
/// a transposition on the wire shows.
fn transform() -> Transform3 {
    Transform3::from_cols_array(&[
        0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 10.0, -20.0, 30.5,
    ])
}

/// One value of every node kind, holding every variant of every enum a node
/// reaches, and the record of which variants it wrote.
fn fixture() -> (GeometryGraph, Coverage) {
    let mut coverage = Coverage::default();
    let mut builder = GeometryGraphBuilder::new();
    let mut push = |node: GeometryNode| builder.push(node).expect("fixture node is valid");
    let mut roots = Vec::new();

    let point2 = push(GeometryNode::Point2(Vec2::new(-0.0, 1.5)));
    let point3 = push(GeometryNode::Point3(Vec3::new(0.1, -0.0, 1e-300)));
    push(GeometryNode::Vector2(Vec2::new(0.6, 0.8)));
    push(GeometryNode::Vector3(Vec3::new(0.0, 0.6, -0.8)));
    push(GeometryNode::Frame2(frame2()));
    push(GeometryNode::Frame3(frame3()));
    push(GeometryNode::Transform(transform()));
    push(GeometryNode::PointList2(vec![
        Vec2::ZERO,
        Vec2::new(1.0, 2.0),
    ]));
    push(GeometryNode::PointList3(vec![
        Vec3::ZERO,
        Vec3::new(1.0, 2.0, 3.0),
    ]));
    push(GeometryNode::BoundingBox(Aabb {
        min: Vec3::new(-1.0, -2.0, -3.0),
        max: Vec3::new(1.0, 2.0, 3.0),
    }));

    let curves2: Vec<NodeId> = curve2s(&mut coverage)
        .into_iter()
        .map(|curve| push(GeometryNode::Curve2(curve)))
        .collect();
    let curves3: Vec<NodeId> = curve3s(&mut coverage)
        .into_iter()
        .map(|curve| push(GeometryNode::Curve3(curve)))
        .collect();
    let surfaces: Vec<NodeId> = surfaces(&mut coverage)
        .into_iter()
        .map(|surface| push(GeometryNode::Surface(surface)))
        .collect();
    let profiles: Vec<NodeId> = profiles(&mut coverage)
        .into_iter()
        .map(|profile| push(GeometryNode::Profile(profile)))
        .collect();
    let solids: Vec<NodeId> = primitives(&mut coverage)
        .into_iter()
        .map(|primitive| push(GeometryNode::Primitive(primitive)))
        .collect();
    let half_space = push(GeometryNode::HalfSpace(HalfSpace {
        boundary: Plane3 {
            origin: Vec3::ZERO,
            normal: Vec3::Z,
        },
        agreement: false,
    }));

    let (line_2d, line_3d, plane) = (curves2[0], curves3[0], surfaces[0]);
    let open_path = push(GeometryNode::Curve2(open_polyline2()));
    let closed_path = push(GeometryNode::Curve2(closed_polyline2()));
    let open_profile = push(GeometryNode::OpenProfile(OpenProfile::new(open_path)));
    let rectangle = profiles[0];

    // Curve relations.
    let transitions = [
        Transition::Discontinuous,
        Transition::Continuous,
        Transition::ContinuousSameGradient,
        Transition::ContinuousSameGradientSameCurvature,
    ];
    coverage.all("Transition", &transitions);
    let selectors = [
        TrimSelector::Parameter(0.0),
        TrimSelector::Point2(Vec2::new(0.0, 0.0)),
        TrimSelector::Point3(Vec3::new(2.0, 0.0, 0.0)),
        TrimSelector::ArcLength(1.5),
    ];
    coverage.all("TrimSelector", &selectors);
    let preferences = [
        TrimmingPreference::Parameter,
        TrimmingPreference::Cartesian,
        TrimmingPreference::Unspecified,
    ];
    coverage.all("TrimmingPreference", &preferences);
    let masters = [
        MasterRepresentation::Curve3d,
        MasterRepresentation::ParameterCurveS1,
        MasterRepresentation::ParameterCurveS2,
        MasterRepresentation::Both,
        MasterRepresentation::Unspecified,
    ];
    coverage.all("MasterRepresentation", &masters);
    let frames = [StationFrame::Section, StationFrame::Plan];
    coverage.all("StationFrame", &frames);

    let mut relations = vec![CurveRelation::Composite {
        segments: transitions
            .iter()
            .enumerate()
            .map(|(index, &transition)| CurveSegment {
                curve: curves2[index],
                same_sense: index % 2 == 0,
                transition,
            })
            .collect(),
    }];
    relations.extend(
        preferences
            .iter()
            .map(|&preference| CurveRelation::Trimmed {
                basis: line_3d,
                start: selectors[..2].to_vec(),
                end: selectors[2..].to_vec(),
                sense_agreement: preference != TrimmingPreference::Cartesian,
                preference,
            }),
    );
    relations.push(CurveRelation::Offset {
        basis: line_2d,
        distance: 0.5,
        reference_direction: None,
    });
    relations.push(CurveRelation::Offset {
        basis: line_3d,
        distance: -0.5,
        reference_direction: Some(Vec3::Z),
    });
    relations.extend(masters.iter().map(|&master| CurveRelation::SurfaceCurve {
        curve_3d: line_3d,
        sides: SurfaceSides::two(plane, line_2d, surfaces[1], curves2[1]),
        master,
    }));
    relations.push(CurveRelation::SurfaceCurve {
        curve_3d: line_3d,
        sides: SurfaceSides::one(plane, line_2d),
        master: MasterRepresentation::Curve3d,
    });
    relations.push(CurveRelation::ParameterCurve {
        basis_surface: plane,
        reference_curve: line_2d,
    });
    relations.extend(frames.iter().map(|&frame| CurveRelation::OffsetByStations {
        basis: line_3d,
        stations: vec![
            Station::at(0.0),
            Station::new(2.0, StationOffsets::new(0.5, -0.25, 0.0)),
        ],
        frame,
    }));
    // Format 1.2: an offset whose stations carry sides and windows (#301),
    // and one whose stations write neither.
    relations.extend(
        frames
            .iter()
            .map(|&frame| CurveRelation::OffsetsAtStations {
                basis: line_3d,
                stations: match frame {
                    StationFrame::Plan => vec![
                        RunStation::new(Station::at(0.0)),
                        RunStation::new(Station::new(2.0, StationOffsets::new(0.5, -0.25, 0.0)))
                            .with_seam_side(SeamSide::Incoming)
                            .with_seam_window(0.001),
                    ],
                    _ => vec![
                        RunStation::new(Station::at(0.5)),
                        RunStation::new(Station::at(1.5)),
                    ],
                },
                frame,
            }),
    );
    // Format 1.3: an offset whose first and last stations stand at its
    // basis's start and end, measured by the kernel (#303).
    let ends = vec![
        RunStation::at_basis_start(StationOffsets::new(0.5, 0.0, 0.0)),
        RunStation::new(Station::new(1.0, StationOffsets::new(0.5, 0.0, 0.0))),
        RunStation::at_basis_end(StationOffsets::new(0.75, 0.0, 0.0)).with_seam_window(0.001),
    ];
    coverage.all(
        "BasisEnd",
        &ends
            .iter()
            .filter_map(|station| station.basis_end)
            .collect::<Vec<_>>(),
    );
    relations.push(CurveRelation::OffsetsAtStations {
        basis: line_3d,
        stations: ends,
        frame: StationFrame::Section,
    });
    coverage.all("CurveRelation", &relations);
    let relations: Vec<NodeId> = relations
        .into_iter()
        .map(|relation| push(GeometryNode::CurveRelation(relation)))
        .collect();

    push(GeometryNode::PointOnCurve(PointOnCurve {
        curve: relations[0],
        parameter: 0.25,
    }));
    push(GeometryNode::PointOnSurface(PointOnSurface {
        surface: plane,
        u: 0.5,
        v: -0.5,
    }));

    // Surface relations.
    let orientation = StationOrientation::new(Some(Vec3::new(0.0, 0.1, 1.0)), None);
    let surface_relations = vec![
        SurfaceRelation::CurveBounded {
            basis: plane,
            boundaries: vec![closed_path, line_3d],
            implicit_outer: true,
        },
        SurfaceRelation::RectangularTrimmed {
            basis: surfaces[1],
            u: (0.0, 1.0),
            v: (-1.0, 1.0),
            u_sense: true,
            v_sense: false,
        },
        SurfaceRelation::Offset {
            basis: plane,
            distance: 0.25,
            self_intersect: Some(false),
        },
        SurfaceRelation::LinearExtrusion {
            swept_curve: line_2d,
            direction: Vec3::Z,
        },
        SurfaceRelation::Revolution {
            swept_curve: line_3d,
            axis_origin: Vec3::ZERO,
            axis_direction: Vec3::Z,
        },
        SurfaceRelation::SectionedSurface {
            directrix: line_3d,
            sections: vec![
                StationedOpenSection {
                    profile: open_profile,
                    tags: vec!["a".into(), "b".into(), "c".into()],
                    station: Station::at(0.0),
                },
                StationedOpenSection {
                    profile: open_profile,
                    tags: vec!["c".into(), "b".into(), "a".into()],
                    station: Station::at(3.0),
                },
            ],
            frame: StationFrame::Plan,
        },
        SurfaceRelation::OpenSectionsAtStations {
            directrix: line_3d,
            sections: vec![
                SectionAtStation::new(open_profile, Station::at(0.0)),
                SectionAtStation::new(open_profile, Station::at(4.0)).with_orientation(orientation),
                // Format 1.2: a window alone (#301).
                SectionAtStation::new(open_profile, Station::at(5.0)).with_seam_window(0.002),
            ],
            frame: StationFrame::Section,
        },
    ];
    coverage.all("SurfaceRelation", &surface_relations);
    for relation in surface_relations {
        push(GeometryNode::SurfaceRelation(relation));
    }

    // Solid operations.
    let operators = [
        BooleanOperator::Union,
        BooleanOperator::Intersection,
        BooleanOperator::Difference,
        BooleanOperator::SymmetricDifference,
    ];
    coverage.all("BooleanOperator", &operators);
    let mut operations = vec![
        SolidOperation::Extrusion {
            profile: rectangle,
            direction: Vec3::Z,
            depth: 3.0,
        },
        SolidOperation::TaperedExtrusion {
            start_profile: rectangle,
            end_profile: profiles[1],
            direction: Vec3::Z,
            depth: 2.0,
        },
        SolidOperation::Revolution {
            profile: rectangle,
            axis_origin: Vec3::new(5.0, 0.0, 0.0),
            axis_direction: Vec3::Y,
            angle: core::f64::consts::FRAC_PI_2,
        },
        SolidOperation::TaperedRevolution {
            start_profile: rectangle,
            end_profile: profiles[2],
            axis_origin: Vec3::ZERO,
            axis_direction: Vec3::Y,
            angle: 1.0,
        },
        SolidOperation::SweptDisk {
            directrix: line_3d,
            radius: 0.05,
            inner_radius: Some(0.04),
            parameter_range: Some((0.0, 2.0)),
            fillet_radius: Some(0.09),
        },
        SolidOperation::SweptDisk {
            directrix: curves3[3],
            radius: 0.05,
            inner_radius: None,
            parameter_range: None,
            fillet_radius: None,
        },
        SolidOperation::FixedReferenceSweep {
            profile: rectangle,
            directrix: line_3d,
            reference_direction: Vec3::Z,
            parameter_range: None,
        },
        SolidOperation::SurfaceCurveSweep {
            profile: rectangle,
            directrix: line_3d,
            reference_surface: plane,
            parameter_range: Some((0.5, 1.5)),
        },
        SolidOperation::SectionedSpine {
            spine: line_3d,
            sections: vec![
                Section {
                    profile: rectangle,
                    placement: Transform3::IDENTITY,
                },
                Section {
                    profile: rectangle,
                    placement: transform(),
                },
            ],
        },
        SolidOperation::BoundedHalfSpace {
            half_space,
            boundary: closed_path,
            placement: transform(),
        },
        SolidOperation::BoundedHalfSpace {
            half_space,
            boundary: profiles[3],
            placement: Transform3::IDENTITY,
        },
        SolidOperation::StationedSpine {
            directrix: line_3d,
            sections: vec![
                StationedSection {
                    profile: rectangle,
                    station: Station::at(0.0),
                },
                StationedSection {
                    profile: rectangle,
                    station: Station::new(5.0, StationOffsets::new(0.0, 1.0, 0.0)),
                },
            ],
            frame: StationFrame::Plan,
        },
        SolidOperation::SectionsAtStations {
            directrix: line_3d,
            sections: vec![
                SectionAtStation::new(rectangle, Station::at(1.0)).with_tags(["p", "q"]),
                SectionAtStation::new(rectangle, Station::at(2.0))
                    .with_tags(["q", "p"])
                    .with_orientation(orientation),
                // Format 1.2: a side alone, and a side with a window (#301).
                SectionAtStation::new(rectangle, Station::at(3.0))
                    .with_tags(["p", "q"])
                    .with_seam_side(SeamSide::Incoming),
                SectionAtStation::new(rectangle, Station::at(4.0))
                    .with_tags(["q", "p"])
                    .with_seam_side(SeamSide::Incoming)
                    .with_seam_window(0.000_25),
            ],
            frame: StationFrame::Section,
        },
    ];
    operations.extend(operators.iter().map(|&operator| SolidOperation::Boolean {
        left: solids[0],
        right: solids[1],
        operator,
    }));
    coverage.all("SolidOperation", &operations);
    for operation in operations {
        roots.push(push(GeometryNode::SolidOperation(operation)));
    }

    // Topology and meshes.
    let mut brep = BRep::<NodeId>::default();
    let start = brep.add_vertex(Vertex {
        position: Vec3::ZERO,
    });
    let end = brep.add_vertex(Vertex { position: Vec3::X });
    let edge = brep.add_edge(Edge {
        start,
        end,
        curve: Some(line_3d),
    });
    brep.add_edge(Edge {
        start: end,
        end: start,
        curve: None,
    });
    let orientations = [Orientation::Forward, Orientation::Reversed];
    coverage.all("Orientation", &orientations);
    let wire = brep.add_loop(Loop {
        edges: vec![
            EdgeUse {
                edge,
                orientation: Orientation::Forward,
                pcurve: Some(line_2d),
            },
            EdgeUse {
                edge,
                orientation: Orientation::Reversed,
                pcurve: None,
            },
        ],
    });
    let face = brep.add_face(Face {
        surface: Some(plane),
        bounds: vec![FaceBound {
            loop_id: wire,
            orientation: Orientation::Forward,
            outer: true,
        }],
        orientation: Orientation::Reversed,
    });
    let shell = brep.add_shell(Shell {
        faces: vec![(face, Orientation::Forward)],
        closed: false,
    });
    let void = brep.add_shell(Shell {
        faces: Vec::new(),
        closed: true,
    });
    brep.add_solid(Solid {
        outer: shell,
        voids: vec![void],
    });
    roots.push(push(GeometryNode::BRep(brep)));

    roots.push(push(GeometryNode::PolygonMesh(PolygonMesh {
        positions: vec![Vec3::ZERO, Vec3::X, Vec3::Y, Vec3::new(1.0, 1.0, 0.0)],
        faces: vec![PolygonFace {
            outer: vec![0, 1, 3, 2],
            holes: vec![vec![]],
        }],
    })));
    let blends = [Blend::Linear, Blend::Nearest, Blend::None];
    coverage.all("Blend", &blends);
    let mut attributes: Vec<AttributeChannel> = blends
        .iter()
        .map(|&blend| AttributeChannel::new(format!("{blend:?}"), vec![0.0, 1.0, 2.0], 1, blend))
        .collect();
    attributes[0].corner_indices = Some(vec![0, 1, 2]);
    roots.push(push(GeometryNode::TriMesh(TriMesh {
        positions: vec![Vec3::ZERO, Vec3::X, Vec3::Y],
        indices: vec![0, 1, 2],
        normals: Some(NormalAttribute {
            values: vec![Vec3::Z],
            indices: Some(vec![0, 0, 0]),
        }),
        attributes,
    })));

    // Instancing and stations.
    let instance = push(GeometryNode::Instance(Instance {
        source: solids[0],
        transform: transform(),
    }));
    let sides = [SeamSide::Outgoing, SeamSide::Incoming];
    coverage.all("SeamSide", &sides);
    let mut station = CurveStation::new(
        line_3d,
        Station::new(1.0, StationOffsets::new(0.5, 0.25, -0.125)),
    );
    station.frame = StationFrame::Plan;
    push(GeometryNode::CurveStation(station));
    let oriented = OrientedCurveStation::new(
        station,
        StationOrientation::new(Some(Vec3::Z), Some(Vec3::new(1.0, 0.5, 0.0))),
    );
    for side in sides {
        push(GeometryNode::OrientedCurveStation(
            oriented.with_seam_side(side),
        ));
    }
    // Format 1.1: a seam-snapping window (#294).
    push(GeometryNode::OrientedCurveStation(
        oriented
            .with_seam_side(SeamSide::Incoming)
            .with_seam_window(0.001),
    ));
    roots.push(push(GeometryNode::InstanceAtStation(
        InstanceAtStation::new(
            line_2d,
            CurveStation::new(curves3[3], Station::at(0.5))
                .with_seam_side(SeamSide::Incoming)
                .with_seam_window(0.000_5),
        ),
    )));
    roots.push(push(GeometryNode::Collection(vec![
        point2, point3, instance,
    ])));

    let graph = builder.finish(roots).expect("fixture roots are valid");
    for (_, node) in graph.iter() {
        coverage.one("GeometryNode", node);
    }
    (graph, coverage)
}

/// Every enum a graph node reaches, with the variants its deserialiser
/// knows.
fn every_enum() -> Vec<(&'static str, BTreeSet<String>)> {
    vec![
        ("GeometryNode", variants_of::<GeometryNode>()),
        ("CurveRelation", variants_of::<CurveRelation>()),
        ("TrimSelector", variants_of::<TrimSelector>()),
        ("TrimmingPreference", variants_of::<TrimmingPreference>()),
        ("Transition", variants_of::<Transition>()),
        (
            "MasterRepresentation",
            variants_of::<MasterRepresentation>(),
        ),
        ("SolidOperation", variants_of::<SolidOperation>()),
        ("SurfaceRelation", variants_of::<SurfaceRelation>()),
        ("StationFrame", variants_of::<StationFrame>()),
        ("SeamSide", variants_of::<SeamSide>()),
        ("BasisEnd", variants_of::<BasisEnd>()),
        ("BooleanOperator", variants_of::<BooleanOperator>()),
        ("Curve2", variants_of::<Curve2>()),
        ("Curve3", variants_of::<Curve3>()),
        ("CurvatureLaw", variants_of::<CurvatureLaw>()),
        ("ChainPiece2", variants_of::<ChainPiece2>()),
        ("ElevationLaw", variants_of::<ElevationLaw>()),
        ("CantForm", variants_of::<CantForm>()),
        ("RailSide", variants_of::<RailSide>()),
        ("BankConvention", variants_of::<BankConvention>()),
        ("KnotSpec", variants_of::<KnotSpec>()),
        ("Basis", variants_of::<Basis>()),
        ("Field2", variants_of::<Field2>()),
        ("Axis", variants_of::<Axis>()),
        ("Carrier", variants_of::<Carrier>()),
        ("Branch", variants_of::<Branch>()),
        ("Surface", variants_of::<Surface>()),
        ("Profile", variants_of::<Profile>()),
        ("SectionProfile", variants_of::<SectionProfile>()),
        ("Primitive", variants_of::<Primitive>()),
        ("Orientation", variants_of::<Orientation>()),
        ("Blend", variants_of::<Blend>()),
    ]
}

/// The fixture writes every variant of every enum a node reaches. A
/// variant added without a fixture value (and so without a golden payload
/// and a minor version) fails here.
#[test]
fn the_fixture_covers_every_variant_of_every_enum() {
    let (graph, coverage) = fixture();
    let text = graph.to_json().unwrap();
    let enums = every_enum();
    let names: BTreeSet<_> = enums.iter().map(|(name, _)| *name).collect();
    let recorded: BTreeSet<_> = coverage.0.keys().copied().collect();
    assert_eq!(
        recorded, names,
        "the coverage record names exactly the enums listed"
    );
    for (name, variants) in enums {
        assert!(!variants.is_empty(), "{name}: no variants probed");
        assert_eq!(
            coverage.0[name], variants,
            "{name}: the fixture must write every variant"
        );
        for variant in &variants {
            assert!(
                text.contains(&format!("\"{variant}\"")),
                "{name}::{variant} is recorded but not in the payload"
            );
        }
    }
}

// ---------------------------------------------------------------------------
// Round trips
// ---------------------------------------------------------------------------

/// The graph part of a payload, keys sorted, numbers in their exact text:
/// a comparison independent of key order and formatting that still tells
/// `-0.0` from `0.0`.
fn canonical_graph(text: &str) -> String {
    let value: Value = serde_json::from_str(text).unwrap();
    serde_json::to_string(&value["graph"]).unwrap()
}

#[test]
fn every_node_kind_round_trips_through_json_and_cbor() {
    let (graph, _) = fixture();
    let text = graph.to_json().unwrap();

    let from_json = GeometryGraph::from_json(&text).unwrap();
    assert_eq!(from_json.len(), graph.len());
    assert_eq!(
        from_json.to_json().unwrap(),
        text,
        "JSON re-encodes identically"
    );

    let bytes = graph.to_cbor().unwrap();
    let from_cbor = GeometryGraph::from_cbor(&bytes).unwrap();
    assert_eq!(
        from_cbor.to_json().unwrap(),
        text,
        "CBOR reads the same graph"
    );
    assert_eq!(
        from_cbor.to_cbor().unwrap(),
        bytes,
        "CBOR re-encodes identically"
    );

    // The rebuilt graph owns its handles: its roots resolve in it, and in no
    // other graph.
    for &root in from_json.roots() {
        assert!(from_json.get(root).is_some());
        assert!(from_cbor.get(root).is_none());
        assert!(graph.get(root).is_none());
    }
    // Payloads are the nodes in order, kind for kind.
    for ((_, before), (_, after)) in graph.iter().zip(from_cbor.iter()) {
        assert_eq!(
            core::mem::discriminant(before),
            core::mem::discriminant(after)
        );
    }
}

#[test]
fn reals_round_trip_bit_exactly() {
    let values = [
        -0.0,
        0.0,
        0.1,
        1.0 / 3.0,
        core::f64::consts::PI,
        f64::MIN_POSITIVE,
        5e-324,
        f64::MAX,
        -f64::MAX,
        1e300,
        -1.5e-10,
        2.0_f64.powi(53) + 2.0,
        0.30000000000000004,
    ];
    let mut builder = GeometryGraphBuilder::new();
    let points: Vec<Point3> = values
        .chunks(3)
        .map(|chunk| Vec3::new(chunk[0], chunk[chunk.len() / 2], chunk[chunk.len() - 1]))
        .collect();
    let list = builder
        .push(GeometryNode::PointList3(points.clone()))
        .unwrap();
    let graph = builder.finish(vec![list]).unwrap();

    let bits = |graph: &GeometryGraph| -> Vec<u64> {
        let (_, GeometryNode::PointList3(points)) = graph.iter().next().unwrap() else {
            panic!("not a point list");
        };
        points
            .iter()
            .flat_map(|point| point.to_array())
            .map(f64::to_bits)
            .collect()
    };
    let expected = bits(&graph);
    let json = GeometryGraph::from_json(&graph.to_json().unwrap()).unwrap();
    assert_eq!(bits(&json), expected, "JSON");
    let cbor = GeometryGraph::from_cbor(&graph.to_cbor().unwrap()).unwrap();
    assert_eq!(bits(&cbor), expected, "CBOR");
    assert!(graph.to_json().unwrap().contains("-0.0"));
}

#[test]
fn many_ordinary_reals_round_trip_bit_exactly() {
    // A JSON reader without correctly rounded parsing lands one ulp off for
    // roughly a quarter of ordinary doubles (12.838242182209141 is one), so a
    // handful of simple values cannot pin the ADR 0085 promise. Sweep a
    // deterministic spread of mantissas across magnitudes instead.
    let mut state: u64 = 0x9E37_79B9_7F4A_7C15;
    let mut values = Vec::new();
    for exponent in -12..=7 {
        for _ in 0..300 {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            let mantissa = 1.0 + (state >> 11) as f64 / (1u64 << 53) as f64;
            let sign = if state & 1 == 0 { 1.0 } else { -1.0 };
            values.push(sign * mantissa * 10f64.powi(exponent));
        }
    }
    values.push(12.838_242_182_209_141);
    values.push(10.003_124_511_871_276);
    while values.len() % 3 != 0 {
        values.push(0.0);
    }
    let points: Vec<Point3> = values
        .chunks(3)
        .map(|c| Vec3::new(c[0], c[1], c[2]))
        .collect();
    let mut builder = GeometryGraphBuilder::new();
    let list = builder.push(GeometryNode::PointList3(points)).unwrap();
    let graph = builder.finish(vec![list]).unwrap();
    let bits = |graph: &GeometryGraph| -> Vec<u64> {
        let (_, GeometryNode::PointList3(points)) = graph.iter().next().unwrap() else {
            panic!("not a point list");
        };
        points
            .iter()
            .flat_map(|point| point.to_array())
            .map(f64::to_bits)
            .collect()
    };
    let expected = bits(&graph);
    let json = GeometryGraph::from_json(&graph.to_json().unwrap()).unwrap();
    let off = expected
        .iter()
        .zip(bits(&json))
        .filter(|(a, b)| **a != *b)
        .count();
    assert_eq!(
        off,
        0,
        "{off} of {} reals changed through JSON",
        values.len()
    );
    let cbor = GeometryGraph::from_cbor(&graph.to_cbor().unwrap()).unwrap();
    assert_eq!(bits(&cbor), expected, "CBOR");
}

// ---------------------------------------------------------------------------
// Golden payloads
// ---------------------------------------------------------------------------

fn golden_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/wire")
}

/// The version entry of a JSON payload.
fn version_of(text: &str) -> String {
    let value: Value = serde_json::from_str(text).unwrap();
    value["version"].as_str().unwrap().to_owned()
}

/// The version entry of a CBOR payload.
fn cbor_version_of(bytes: &[u8]) -> String {
    let value: ciborium::Value = ciborium::from_reader(bytes).unwrap();
    let ciborium::Value::Map(entries) = value else {
        panic!("a CBOR payload is a map")
    };
    entries
        .iter()
        .find(|(key, _)| key.as_text() == Some("version"))
        .and_then(|(_, version)| version.as_text())
        .unwrap()
        .to_owned()
}

/// A golden payload stored one node per line, as the writer writes it: the
/// line breaks between nodes are the only bytes the writer does not write.
fn as_written(golden: &str) -> String {
    golden.replace('\n', "")
}

/// A payload with one node per line, as the golden files are stored, under
/// the version the writer labels it with.
fn one_node_per_line(graph: &GeometryGraph) -> String {
    let nodes: Vec<String> = graph
        .iter()
        .map(|(_, node)| serde_json::to_string(node).unwrap())
        .collect();
    let roots: Vec<String> = graph
        .roots()
        .iter()
        .map(|root| root.index().to_string())
        .collect();
    let version = version_of(&graph.to_json().unwrap());
    format!(
        "{{\"format\":\"{FORMAT_NAME}\",\"version\":\"{version}\",\"graph\":{{\"nodes\":[\n{}\n],\n\"roots\":[{}]}}}}\n",
        nodes.join(",\n"),
        roots.join(",")
    )
}

/// Every published format version's golden payloads read, and write back
/// the same graph. This is the stability promise's regression guard: a
/// change that breaks it is a major version.
#[test]
fn every_golden_payload_reads_to_the_graph_it_was_written_from() {
    let mut read = 0;
    for version in std::fs::read_dir(golden_dir()).unwrap() {
        let version = version.unwrap().path();
        for payload in std::fs::read_dir(&version).unwrap() {
            let path = payload.unwrap().path();
            let text = std::fs::read_to_string(&path).unwrap();
            let graph = GeometryGraph::from_json(&text)
                .unwrap_or_else(|error| panic!("{}: {error}", path.display()));
            assert_eq!(
                canonical_graph(&graph.to_json().unwrap()),
                canonical_graph(&text),
                "{} reads back changed",
                path.display()
            );
            let cbor = GeometryGraph::from_cbor(&graph.to_cbor().unwrap()).unwrap();
            assert_eq!(cbor.to_json().unwrap(), graph.to_json().unwrap());
            read += 1;
        }
    }
    assert!(read > 0, "no golden payloads found");
}

/// Today's writer writes the current version's golden payload exactly. Set
/// `AXIOLID_BLESS_WIRE_GOLDEN=1` to write it when a NEW format version is
/// introduced; a published version's payload is never rewritten.
#[test]
fn the_writer_writes_the_current_golden_payload() {
    let (graph, _) = fixture();
    let path = golden_dir()
        .join(format!("v{FORMAT_VERSION}"))
        .join("every-kind.json");
    let written = one_node_per_line(&graph);
    if std::env::var_os("AXIOLID_BLESS_WIRE_GOLDEN").is_some() {
        std::fs::create_dir_all(path.parent().unwrap()).unwrap();
        std::fs::write(&path, &written).unwrap();
    }
    let golden = std::fs::read_to_string(&path).unwrap();
    assert_eq!(
        written,
        golden,
        "the writer's output for {}",
        path.display()
    );
    assert_eq!(graph.to_json().unwrap(), as_written(&golden));
    assert_eq!(
        cbor_version_of(&graph.to_cbor().unwrap()),
        FORMAT_VERSION.to_string()
    );
}

/// Every published version's every-kind payload is written back byte for
/// byte, under its own version (#297): a graph that uses nothing a later
/// minor added is written exactly as that version's writer wrote it, so
/// its readers keep reading it. The 1.0 payload is what 0.3.8 wrote.
#[test]
fn every_golden_payload_is_written_back_identically_under_its_version() {
    let mut written = Vec::new();
    for version in std::fs::read_dir(golden_dir()).unwrap() {
        let version = version.unwrap().path();
        let path = version.join("every-kind.json");
        let golden = std::fs::read_to_string(&path).unwrap();
        let graph = GeometryGraph::from_json(&golden).unwrap();
        let label = version_of(&golden);
        assert_eq!(
            format!("v{label}"),
            version.file_name().unwrap().to_str().unwrap()
        );
        assert_eq!(
            graph.to_json().unwrap(),
            as_written(&golden),
            "{} is not written back byte for byte",
            path.display()
        );
        assert_eq!(cbor_version_of(&graph.to_cbor().unwrap()), label, "CBOR");
        written.push(label);
    }
    written.sort();
    assert_eq!(written, ["1.0", "1.1", "1.2", "1.3"]);
}

// ---------------------------------------------------------------------------
// Refusals
// ---------------------------------------------------------------------------

fn small_graph() -> GeometryGraph {
    let mut builder = GeometryGraphBuilder::new();
    let sphere = builder
        .push(GeometryNode::Primitive(Primitive::Sphere { radius: 1.0 }))
        .unwrap();
    let instance = builder
        .push(GeometryNode::Instance(Instance {
            source: sphere,
            transform: Transform3::IDENTITY,
        }))
        .unwrap();
    builder.finish(vec![instance]).unwrap()
}

/// The version the small graph is written under: it has no window, so the
/// lowest version that represents it (#297).
const ONE_ZERO: FormatVersion = FormatVersion::new(1, 0);

fn small_payload() -> Value {
    serde_json::from_str(&small_graph().to_json().unwrap()).unwrap()
}

fn read_json(value: &Value) -> Result<GeometryGraph, WireError> {
    GeometryGraph::from_json(&value.to_string())
}

fn read_cbor(value: &Value) -> Result<GeometryGraph, WireError> {
    let mut bytes = Vec::new();
    ciborium::into_writer(value, &mut bytes).unwrap();
    GeometryGraph::from_cbor(&bytes)
}

/// Both encodings refuse `value` with `expected`.
fn refused(value: &Value, expected: &WireError) {
    assert_eq!(read_json(value).unwrap_err(), *expected, "JSON: {value}");
    assert_eq!(read_cbor(value).unwrap_err(), *expected, "CBOR: {value}");
}

#[test]
fn the_small_payload_has_the_documented_shape() {
    assert_eq!(
        small_payload(),
        json!({
            "format": "axiolid-geometry-graph",
            "version": "1.0",
            "graph": {
                "nodes": [
                    {"Primitive": {"Sphere": {"radius": 1.0}}},
                    {"Instance": {
                        "source": 0,
                        "transform": [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0]
                    }}
                ],
                "roots": [1]
            }
        })
    );
    assert!(read_json(&small_payload()).is_ok());
    assert!(read_cbor(&small_payload()).is_ok());
}

#[test]
fn an_unknown_node_kind_is_refused_by_name_and_version() {
    let mut payload = small_payload();
    payload["graph"]["nodes"][0] = json!({"Hyperboloid": {"radius": 1.0}});
    refused(
        &payload,
        &WireError::UnknownKind {
            kind: "Hyperboloid".into(),
            version: ONE_ZERO,
        },
    );
}

#[test]
fn an_unknown_variant_is_refused_by_name() {
    let mut payload = small_payload();
    payload["graph"]["nodes"][0] = json!({"Primitive": {"Ellipsoid": {"radius": 1.0}}});
    refused(
        &payload,
        &WireError::UnknownKind {
            kind: "Ellipsoid".into(),
            version: ONE_ZERO,
        },
    );
    // A unit variant too.
    let (graph, _) = fixture();
    let text = graph.to_json().unwrap();
    assert!(text.contains("\"SymmetricDifference\""));
    let edited = text.replacen("\"SymmetricDifference\"", "\"Xor\"", 1);
    assert_eq!(
        GeometryGraph::from_json(&edited).unwrap_err(),
        WireError::UnknownKind {
            kind: "Xor".into(),
            version: FORMAT_VERSION,
        }
    );
}

#[test]
fn an_unknown_field_is_refused_by_name() {
    let mut payload = small_payload();
    payload["graph"]["nodes"][0]["Primitive"]["Sphere"]["colour"] = json!("red");
    let expected = WireError::UnknownField {
        field: "colour".into(),
        version: ONE_ZERO,
    };
    refused(&payload, &expected);
    let mut payload = small_payload();
    payload["graph"]["colour"] = json!(1);
    refused(&payload, &expected);
    let mut payload = small_payload();
    payload["colour"] = json!(1);
    refused(&payload, &expected);
}

#[test]
fn a_version_this_reader_does_not_read_is_refused_before_the_graph() {
    for found in ["2.0", "1.4", "0.9", "1", "one", "1.0.0", ""] {
        let mut payload = small_payload();
        payload["version"] = json!(found);
        // Even with content this reader does not know: the version is
        // checked first, so a newer payload is never misreported.
        payload["graph"]["nodes"][0] = json!({"Hyperboloid": {}});
        refused(
            &payload,
            &WireError::UnsupportedVersion {
                found: found.into(),
                supported: FORMAT_VERSION,
            },
        );
    }
}

#[test]
fn a_missing_version_or_another_format_is_refused() {
    let mut payload = small_payload();
    payload.as_object_mut().unwrap().remove("version");
    refused(&payload, &WireError::MissingVersion);

    let mut payload = small_payload();
    payload.as_object_mut().unwrap().remove("format");
    refused(&payload, &WireError::UnknownFormat { found: None });

    let mut payload = small_payload();
    payload["format"] = json!("geometry");
    refused(
        &payload,
        &WireError::UnknownFormat {
            found: Some("geometry".into()),
        },
    );
}

#[test]
fn a_dangling_or_forward_reference_is_refused_by_validation() {
    for source in [1, 2, 99] {
        let mut payload = small_payload();
        payload["graph"]["nodes"][1]["Instance"]["source"] = json!(source);
        let error = read_json(&payload).unwrap_err();
        assert!(
            matches!(
                error,
                WireError::InvalidGraph {
                    node: Some(1),
                    error: GraphError::NonPriorReference { .. },
                }
            ),
            "{source}: {error:?}"
        );
        // The same refusal; the handles in it carry each read's own brand.
        assert_eq!(
            read_cbor(&payload).unwrap_err().to_string(),
            error.to_string()
        );
    }

    let mut payload = small_payload();
    payload["graph"]["roots"] = json!([2]);
    assert!(matches!(
        read_json(&payload).unwrap_err(),
        WireError::InvalidGraph {
            node: None,
            error: GraphError::UnknownRoot { node_count: 2, .. },
        }
    ));
}

#[test]
fn a_reference_to_the_wrong_family_is_refused_by_validation() {
    let mut payload = small_payload();
    payload["graph"]["nodes"][1] = json!({"SolidOperation": {"Extrusion": {
        "profile": 0, "direction": [0.0, 0.0, 1.0], "depth": 1.0
    }}});
    assert!(matches!(
        read_json(&payload).unwrap_err(),
        WireError::InvalidGraph {
            node: Some(1),
            error: GraphError::InvalidReferenceType {
                expected: "profile",
                actual: "primitive",
                ..
            },
        }
    ));
}

#[test]
fn an_invalid_station_is_refused_by_validation() {
    let mut payload = small_payload();
    payload["graph"]["nodes"][0] = json!({"Curve3": {"Line": {
        "origin": [0.0, 0.0, 0.0], "direction": [1.0, 0.0, 0.0]
    }}});
    payload["graph"]["nodes"][1] = json!({"CurveStation": {
        "basis": 0,
        "station": {"distance": -1.0, "offsets": {"lateral": 0.0, "vertical": 0.0, "longitudinal": 0.0}},
        "frame": "Section"
    }});
    assert!(matches!(
        read_json(&payload).unwrap_err(),
        WireError::InvalidGraph {
            node: Some(1),
            error: GraphError::InvalidStation { .. },
        }
    ));
}

#[test]
fn a_non_finite_number_is_refused_on_read_and_write() {
    // CBOR can carry a NaN or an infinity; the reader refuses it with its
    // path.
    fn replace(value: &mut ciborium::Value, marker: f64, bad: f64) {
        match value {
            ciborium::Value::Float(number) if *number == marker => *number = bad,
            ciborium::Value::Array(items) => {
                items.iter_mut().for_each(|item| replace(item, marker, bad));
            }
            ciborium::Value::Map(entries) => entries
                .iter_mut()
                .for_each(|(_, item)| replace(item, marker, bad)),
            _ => {}
        }
    }
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let mut payload = small_payload();
        payload["graph"]["nodes"][0]["Primitive"]["Sphere"]["radius"] = json!(12345.5);
        let mut value = ciborium::Value::serialized(&payload).unwrap();
        replace(&mut value, 12345.5, bad);
        let mut bytes = Vec::new();
        ciborium::into_writer(&value, &mut bytes).unwrap();
        assert_eq!(
            GeometryGraph::from_cbor(&bytes).unwrap_err(),
            WireError::NonFinite {
                path: "graph.nodes[0].Primitive.Sphere.radius".into()
            }
        );
    }

    // JSON has no literal for one; a literal past the double range is one.
    let text = small_payload()
        .to_string()
        .replace("\"radius\":1.0", "\"radius\":1e999");
    assert!(matches!(
        GeometryGraph::from_json(&text).unwrap_err(),
        WireError::NonFinite { .. }
    ));

    // A writer refuses rather than write `null`.
    let mut builder = GeometryGraphBuilder::new();
    let sphere = builder
        .push(GeometryNode::Primitive(Primitive::Sphere {
            radius: f64::NAN,
        }))
        .unwrap();
    let point = builder
        .push(GeometryNode::Point3(Vec3::new(0.0, f64::INFINITY, 0.0)))
        .unwrap();
    let graph = builder.finish(vec![sphere, point]).unwrap();
    let expected = WireError::NonFinite {
        path: "graph.nodes[0].Primitive.Sphere.radius".into(),
    };
    assert_eq!(graph.to_json().unwrap_err(), expected);
    assert_eq!(graph.to_cbor().unwrap_err(), expected);
}

#[test]
fn malformed_payloads_are_refused() {
    let bytes = small_graph().to_cbor().unwrap();
    for broken in [
        &bytes[..bytes.len() - 1],
        &bytes[..1],
        &[][..],
        &[0xff, 0x00, 0x13][..],
    ] {
        assert!(
            matches!(
                GeometryGraph::from_cbor(broken).unwrap_err(),
                WireError::Malformed { .. }
            ),
            "{broken:02x?}"
        );
    }
    let mut trailing = bytes.clone();
    trailing.push(0);
    assert!(matches!(
        GeometryGraph::from_cbor(&trailing).unwrap_err(),
        WireError::Malformed { .. }
    ));

    for text in ["", "[]", "{", "null", "{\"format\":3}"] {
        assert!(
            matches!(
                GeometryGraph::from_json(text).unwrap_err(),
                WireError::Malformed { .. }
            ),
            "{text:?}"
        );
    }

    // A required field missing, and a wrong type.
    let mut payload = small_payload();
    payload["graph"]["nodes"][1]["Instance"]
        .as_object_mut()
        .unwrap()
        .remove("transform");
    assert!(matches!(
        read_json(&payload).unwrap_err(),
        WireError::Malformed { .. }
    ));
    let mut payload = small_payload();
    payload["graph"]["nodes"][0]["Primitive"]["Sphere"]["radius"] = json!("1.0");
    assert!(matches!(
        read_cbor(&payload).unwrap_err(),
        WireError::Malformed { .. }
    ));
}

/// The small payload as CBOR with the sphere's radius written as `radius`,
/// any CBOR item.
fn cbor_with_radius(radius: ciborium::Value) -> Vec<u8> {
    fn replace(value: &mut ciborium::Value, radius: &ciborium::Value) {
        match value {
            ciborium::Value::Float(number) if *number == 12345.5 => *value = radius.clone(),
            ciborium::Value::Array(items) => {
                items.iter_mut().for_each(|item| replace(item, radius));
            }
            ciborium::Value::Map(entries) => entries
                .iter_mut()
                .for_each(|(_, item)| replace(item, radius)),
            _ => {}
        }
    }
    let mut payload = small_payload();
    payload["graph"]["nodes"][0]["Primitive"]["Sphere"]["radius"] = json!(12345.5);
    let mut value = ciborium::Value::serialized(&payload).unwrap();
    replace(&mut value, &radius);
    let mut bytes = Vec::new();
    ciborium::into_writer(&value, &mut bytes).unwrap();
    bytes
}

fn sphere_radius(graph: &GeometryGraph) -> f64 {
    let (_, GeometryNode::Primitive(Primitive::Sphere { radius })) = graph.iter().next().unwrap()
    else {
        panic!("not a sphere");
    };
    *radius
}

/// An integer where a real is expected reads as that real when it is
/// exactly a double, in both encodings: host encoders write whole-number
/// doubles as integers.
#[test]
fn an_exact_integer_reads_as_a_real() {
    let mut payload = small_payload();
    payload["graph"]["nodes"][0]["Primitive"]["Sphere"]["radius"] = json!(1);
    assert_eq!(sphere_radius(&read_json(&payload).unwrap()), 1.0);

    let two_53 = 1_i128 << 53;
    for (integer, expected) in [
        (1_i128, 1.0),
        (-3, -3.0),
        (two_53, 9_007_199_254_740_992.0),
        (-two_53, -9_007_199_254_740_992.0),
        // Past 2^53 an integer still reads when it is a double.
        (1 << 60, 1_152_921_504_606_846_976.0),
        // CBOR's most negative integer, -2^64, is a double too.
        (-i128::from(u64::MAX) - 1, -18_446_744_073_709_551_616.0),
    ] {
        let bytes = cbor_with_radius(ciborium::Value::Integer(integer.try_into().unwrap()));
        let radius = sphere_radius(&GeometryGraph::from_cbor(&bytes).unwrap());
        assert_eq!(radius.to_bits(), f64::to_bits(expected), "{integer}");
    }
}

/// Integer zero has no sign: CBOR has no integer `-0`, so a host that
/// writes `-0.0` as an integer zero reads back `+0.0`. A float `-0.0` keeps
/// its sign.
#[test]
fn integer_zero_reads_as_positive_zero() {
    let bytes = cbor_with_radius(ciborium::Value::Integer(0.into()));
    let radius = sphere_radius(&GeometryGraph::from_cbor(&bytes).unwrap());
    assert_eq!(radius.to_bits(), 0.0_f64.to_bits());
    let bytes = cbor_with_radius(ciborium::Value::Float(-0.0));
    let radius = sphere_radius(&GeometryGraph::from_cbor(&bytes).unwrap());
    assert_eq!(radius.to_bits(), (-0.0_f64).to_bits());
}

/// An integer that is not exactly a double is refused by name and path,
/// never rounded, wherever it stands.
#[test]
fn an_inexact_integer_is_refused_with_its_path() {
    let two_53 = 1_i128 << 53;
    for integer in [
        two_53 + 1,
        -(two_53 + 1),
        i128::from(u64::MAX),
        -i128::from(u64::MAX),
    ] {
        let bytes = cbor_with_radius(ciborium::Value::Integer(integer.try_into().unwrap()));
        let error = GeometryGraph::from_cbor(&bytes).unwrap_err();
        let WireError::Malformed { detail } = &error else {
            panic!("{integer}: {error:?}");
        };
        assert!(
            detail.contains("not exactly a double")
                && detail.contains("graph.nodes[0].Primitive.Sphere.radius"),
            "{integer}: {detail}"
        );
    }
    // An index is no exception.
    let mut value = ciborium::Value::serialized(&small_payload()).unwrap();
    let ciborium::Value::Map(entries) = &mut value else {
        panic!("a payload is a map")
    };
    for (_, item) in entries.iter_mut() {
        if let ciborium::Value::Map(graph) = item {
            for (key, roots) in graph.iter_mut() {
                if key.as_text() == Some("roots") {
                    *roots = ciborium::Value::Array(vec![ciborium::Value::Integer(
                        u64::try_from(two_53 + 1).unwrap().into(),
                    )]);
                }
            }
        }
    }
    let mut bytes = Vec::new();
    ciborium::into_writer(&value, &mut bytes).unwrap();
    assert!(matches!(
        GeometryGraph::from_cbor(&bytes).unwrap_err(),
        WireError::Malformed { detail } if detail.contains("graph.roots[0]")
    ));
}

/// Only the self-describe tag is unwrapped; byte strings and other tags
/// are not part of the format.
#[test]
fn cbor_items_outside_the_data_model_are_refused() {
    let tagged = cbor_with_radius(ciborium::Value::Tag(
        55799,
        Box::new(ciborium::Value::Float(2.0)),
    ));
    assert_eq!(
        sphere_radius(&GeometryGraph::from_cbor(&tagged).unwrap()),
        2.0
    );
    for item in [
        ciborium::Value::Tag(1, Box::new(ciborium::Value::Float(2.0))),
        ciborium::Value::Bytes(vec![1, 2]),
    ] {
        assert!(matches!(
            GeometryGraph::from_cbor(&cbor_with_radius(item)).unwrap_err(),
            WireError::Malformed { .. }
        ));
    }
}

/// An absent optional field reads as absent: what lets a minor version add
/// one while older payloads stay readable.
#[test]
fn a_missing_optional_field_reads_as_absent() {
    let mut payload = small_payload();
    payload["graph"]["nodes"][0] = json!({"Profile": {"Rectangle": {"x": 1.0, "y": 2.0}}});
    payload["graph"]["nodes"][1] = json!({"Curve2": {"Line": {
        "origin": [0.0, 0.0], "direction": [1.0, 0.0]
    }}});
    let graph = read_json(&payload).unwrap();
    let (_, GeometryNode::Profile(Profile::Rectangle(rectangle))) = graph.iter().next().unwrap()
    else {
        panic!("not a rectangle");
    };
    assert_eq!(rectangle.thickness, None);
    assert_eq!(rectangle.outer_radius, None);
    read_cbor(&payload).unwrap();
}

#[test]
fn a_node_reference_does_not_deserialise_outside_a_graph_payload() {
    assert!(serde_json::from_str::<NodeId>("0").is_err());
    assert!(serde_json::from_str::<GeometryNode>(
        r#"{"Instance": {"source": 0, "transform": [1,0,0,0,1,0,0,0,1,0,0,0]}}"#
    )
    .is_err());
    // A node with no reference needs no graph.
    assert!(serde_json::from_str::<GeometryNode>(r#"{"Point2": [1.0, 2.0]}"#).is_ok());
}

#[test]
fn the_format_version_is_named_in_every_payload() {
    assert_eq!(FORMAT_VERSION, FormatVersion::new(1, 3));
    assert_eq!(FORMAT_VERSION.to_string(), "1.3");
    let payload = small_payload();
    assert_eq!(payload["format"], FORMAT_NAME);
    // A graph without windows is labelled 1.0, the lowest version that
    // represents it (#297).
    assert_eq!(payload["version"], "1.0");
    let value: ciborium::Value =
        ciborium::from_reader(&small_graph().to_cbor().unwrap()[..]).unwrap();
    let ciborium::Value::Map(entries) = value else {
        panic!("a CBOR payload is a map")
    };
    let keys: Vec<_> = entries
        .iter()
        .map(|(key, _)| key.as_text().unwrap().to_owned())
        .collect();
    assert_eq!(keys, ["format", "version", "graph"]);
    // The free functions and the methods are one API.
    assert_eq!(
        wire::to_json(&small_graph()).unwrap(),
        small_graph().to_json().unwrap()
    );
}

// ---------------------------------------------------------------------------
// Format 1.1: the seam-snapping window (#294)
// ---------------------------------------------------------------------------

/// A payload holding one line and one oriented station along it, the
/// station's JSON object given.
fn station_payload(version: &str, station: &Value) -> Value {
    json!({
        "format": FORMAT_NAME,
        "version": version,
        "graph": {
            "nodes": [
                {"Curve3": {"Line": {"origin": [0.0, 0.0, 0.0], "direction": [1.0, 0.0, 0.0]}}},
                {"OrientedCurveStation": station},
            ],
            "roots": [1],
        }
    })
}

fn oriented_station(graph: &GeometryGraph) -> OrientedCurveStation {
    match graph.iter().nth(1) {
        Some((_, GeometryNode::OrientedCurveStation(station))) => *station,
        other => panic!("not an oriented station: {other:?}"),
    }
}

fn station_object(window: Option<Value>) -> Value {
    let mut station = json!({
        "station": {
            "basis": 0,
            "station": {
                "distance": 2.0,
                "offsets": {"lateral": 0.0, "vertical": 0.0, "longitudinal": 0.0}
            },
            "frame": "Section",
        },
        "orientation": {"axis": null, "ref_direction": null},
        "seam": "Incoming",
    });
    if let Some(window) = window {
        station["seam_window"] = window;
    }
    station
}

/// A 1.0 payload has no window and reads as a station without one, which
/// is what it meant; today's reader still reads it.
#[test]
fn a_one_zero_station_reads_with_no_window() {
    for version in ["1.0", "1.1"] {
        let payload = station_payload(version, &station_object(None));
        for graph in [read_json(&payload).unwrap(), read_cbor(&payload).unwrap()] {
            let station = oriented_station(&graph);
            assert_eq!(station.seam_window.to_bits(), 0.0_f64.to_bits());
            assert_eq!(station.seam, SeamSide::Incoming);
        }
    }
}

/// A window is a 1.1 field: a payload claiming 1.0 that holds one is
/// refused by name, as any field its version does not have.
#[test]
fn a_window_in_a_one_zero_payload_is_refused_by_name() {
    let payload = station_payload("1.0", &station_object(Some(json!(0.001))));
    refused(
        &payload,
        &WireError::UnknownField {
            field: "seam_window".into(),
            version: FormatVersion::new(1, 0),
        },
    );
    let graph = read_json(&station_payload("1.1", &station_object(Some(json!(0.001))))).unwrap();
    assert_eq!(oriented_station(&graph).seam_window, 0.001);
}

/// A window round-trips bit-exactly through both encodings, a zero window
/// is not written (so a graph without windows writes the 1.0 graph), and a
/// negative or non-finite one is refused by validation.
#[test]
fn a_window_round_trips_and_a_zero_window_is_not_written() {
    let windows = [0.0, -0.0, 1e-3, 0.1 + 0.2];
    let mut builder = GeometryGraphBuilder::new();
    let line = builder.push(GeometryNode::Curve3(line3())).unwrap();
    let plain = OrientedCurveStation::from(CurveStation::new(line, Station::at(2.0)));
    let mut ids = Vec::new();
    for window in windows {
        ids.push(
            builder
                .push(GeometryNode::OrientedCurveStation(
                    plain.with_seam_window(window),
                ))
                .unwrap(),
        );
    }
    let graph = builder.finish(ids).unwrap();
    let text = graph.to_json().unwrap();
    assert_eq!(text.matches("seam_window").count(), 3, "{text}");
    let want: Vec<u64> = windows.iter().map(|w| w.to_bits()).collect();
    for read in [
        GeometryGraph::from_json(&text).unwrap(),
        GeometryGraph::from_cbor(&graph.to_cbor().unwrap()).unwrap(),
    ] {
        let got: Vec<u64> = read
            .iter()
            .filter_map(|(_, node)| match node {
                GeometryNode::OrientedCurveStation(s) => Some(s.seam_window.to_bits()),
                _ => None,
            })
            .collect();
        assert_eq!(got, want);
    }
    for bad in [-1e-3, f64::NAN, f64::INFINITY] {
        let mut builder = GeometryGraphBuilder::new();
        let line = builder.push(GeometryNode::Curve3(line3())).unwrap();
        let station = OrientedCurveStation::from(CurveStation::new(line, Station::at(2.0)))
            .with_seam_window(bad);
        assert!(matches!(
            builder.push(GeometryNode::OrientedCurveStation(station)),
            Err(GraphError::InvalidStation { .. })
        ));
        assert!(matches!(
            builder.push(GeometryNode::InstanceAtStation(InstanceAtStation::new(
                line, station
            ))),
            Err(GraphError::InvalidStation { .. })
        ));
    }
}

// ---------------------------------------------------------------------------
// The lowest version that represents a graph (#297)
// ---------------------------------------------------------------------------

/// Whether a JSON object is an oriented curve station, recognised by its
/// fields wherever it stands, so that a station a new node kind (or a new
/// field of an old one) comes to carry is found without naming it here.
fn is_oriented_station(value: &Value) -> bool {
    value.as_object().is_some_and(|object| {
        ["station", "orientation", "seam"]
            .iter()
            .all(|field| object.contains_key(*field))
            && object.keys().all(|key| {
                ["station", "orientation", "seam", "seam_window"].contains(&key.as_str())
            })
    })
}

/// The JSON pointer of every oriented station inside `value`, depth first.
fn station_pointers(value: &Value, at: &str, found: &mut Vec<String>) {
    if is_oriented_station(value) {
        found.push(at.to_owned());
    }
    match value {
        Value::Object(object) => {
            for (key, item) in object {
                station_pointers(item, &format!("{at}/{key}"), found);
            }
        }
        Value::Array(items) => {
            for (index, item) in items.iter().enumerate() {
                station_pointers(item, &format!("{at}/{index}"), found);
            }
        }
        _ => {}
    }
}

/// `value` without any `seam_window` entry.
fn without_windows(value: &mut Value) {
    match value {
        Value::Object(object) => {
            object.remove("seam_window");
            object.values_mut().for_each(without_windows);
        }
        Value::Array(items) => items.iter_mut().for_each(without_windows),
        _ => {}
    }
}

/// A graph is labelled with the lowest version that represents it: 1.0
/// unless a station anywhere in it carries a window, then 1.1. Every
/// oriented station in the every-kind fixture, which holds every node kind
/// (see `the_fixture_covers_every_variant_of_every_enum`), is given a
/// window in turn, alone; each one must lift the payload to 1.1, and a 1.0
/// reader must refuse that payload, so the writer never labels a payload
/// with a version whose reader would refuse it. What format 1.2 added is
/// taken out of the fixture first (`without_run_seams`).
#[test]
fn a_window_in_any_station_carrier_and_only_a_window_selects_one_one() {
    let (graph, _) = fixture();
    let mut payload: Value = serde_json::from_str(&graph.to_json().unwrap()).unwrap();
    without_run_seams(&mut payload);
    without_windows(&mut payload);
    payload["version"] = json!("1.0");

    let plain = read_json(&payload).unwrap();
    let text = plain.to_json().unwrap();
    assert_eq!(version_of(&text), "1.0");
    assert!(!text.contains("seam_window"));
    assert_eq!(cbor_version_of(&plain.to_cbor().unwrap()), "1.0");
    let zero_windows = graph.iter().filter(|(_, node)| {
        matches!(
            node,
            GeometryNode::OrientedCurveStation(station) if station.seam_window == 0.0
        )
    });
    assert!(
        zero_windows.count() > 0,
        "the fixture has window-free stations"
    );

    let mut pointers = Vec::new();
    station_pointers(&payload["graph"]["nodes"], "", &mut pointers);
    let carriers: BTreeSet<String> = pointers
        .iter()
        .map(|pointer| {
            let index: usize = pointer.split('/').nth(1).unwrap().parse().unwrap();
            tag(&payload["graph"]["nodes"][index])
        })
        .collect();
    assert_eq!(
        carriers,
        BTreeSet::from(["InstanceAtStation".into(), "OrientedCurveStation".into()]),
        "every node kind that carries an oriented station"
    );

    for pointer in &pointers {
        for window in [json!(0.25), json!(-0.0)] {
            let mut windowed = payload.clone();
            windowed["version"] = json!("1.1");
            windowed["graph"]["nodes"]
                .pointer_mut(pointer)
                .unwrap()
                .as_object_mut()
                .unwrap()
                .insert("seam_window".into(), window.clone());
            let graph = read_json(&windowed).unwrap();
            let text = graph.to_json().unwrap();
            assert_eq!(version_of(&text), "1.1", "a window {window} at {pointer}");
            assert_eq!(
                cbor_version_of(&graph.to_cbor().unwrap()),
                "1.1",
                "CBOR: a window {window} at {pointer}"
            );
            windowed["version"] = json!("1.0");
            refused(
                &windowed,
                &WireError::UnknownField {
                    field: "seam_window".into(),
                    version: FormatVersion::new(1, 0),
                },
            );
        }
    }
}

// ---------------------------------------------------------------------------
// Format 1.2: seam sides and windows on run stations (#301)
// ---------------------------------------------------------------------------

/// `value` with what format 1.3 added taken out (#303): every run station
/// at its basis's start stated at `0`, which it is, and every one at its
/// basis's end stated one past its predecessor. The fixture's 1.2 content
/// stays as it is.
fn without_basis_ends(value: &mut Value) {
    match value {
        Value::Object(object) => {
            if let Some(stations) = object
                .get_mut("OffsetsAtStations")
                .and_then(|offsets| offsets["stations"].as_array_mut())
            {
                let mut previous = 0.0;
                for station in stations {
                    let station = station.as_object_mut().unwrap();
                    if station.remove("basis_end") == Some(json!("End")) {
                        station["station"]["distance"] = json!(previous + 1.0);
                    }
                    previous = station["station"]["distance"].as_f64().unwrap();
                }
            }
            object.values_mut().for_each(without_basis_ends);
        }
        Value::Array(items) => items.iter_mut().for_each(without_basis_ends),
        _ => {}
    }
}

/// `value` with what formats 1.2 and 1.3 added taken out: every
/// `OffsetsAtStations` written as the `OffsetByStations` of its plain
/// stations (its basis's ends stated, [`without_basis_ends`]), and every
/// section at a station without its side and window. The fixture's 1.1
/// content stays as it is.
fn without_run_seams(value: &mut Value) {
    without_basis_ends(value);
    without_side_or_window(value);
}

fn without_side_or_window(value: &mut Value) {
    match value {
        Value::Object(object) => {
            if let Some(mut offsets) = object.remove("OffsetsAtStations") {
                let stations: Vec<Value> = offsets["stations"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|station| station["station"].clone())
                    .collect();
                offsets["stations"] = Value::Array(stations);
                object.insert("OffsetByStations".into(), offsets);
            }
            if object.contains_key("profile") && object.contains_key("orientation") {
                object.remove("seam");
                object.remove("seam_window");
            }
            object.values_mut().for_each(without_side_or_window);
        }
        Value::Array(items) => items.iter_mut().for_each(without_side_or_window),
        _ => {}
    }
}

/// A node's kind as the wire names it, with the variant of a relation or
/// an operation: `SolidOperation::SectionsAtStations`, `CurveStation`.
fn kind_of(node: &Value) -> String {
    let outer = tag(node);
    match node[&outer].as_object() {
        Some(inner) if inner.len() == 1 => {
            let variant = inner.keys().next().unwrap();
            if variant.starts_with(char::is_uppercase) {
                format!("{outer}::{variant}")
            } else {
                outer
            }
        }
        _ => outer,
    }
}

/// The JSON pointer of every object inside `value` that is a plain station
/// (a distance and offsets) or holds one as its `station`, depth first:
/// every place a run station's side and window could be written.
fn plain_station_pointers(value: &Value, at: &str, found: &mut Vec<String>) {
    if let Some(object) = value.as_object() {
        let is_station = object.contains_key("distance") && object.contains_key("offsets");
        let holds_station = object
            .get("station")
            .and_then(Value::as_object)
            .is_some_and(|station| station.contains_key("distance"));
        if is_station || holds_station {
            found.push(at.to_owned());
        }
    }
    match value {
        Value::Object(object) => {
            for (key, item) in object {
                plain_station_pointers(item, &format!("{at}/{key}"), found);
            }
        }
        Value::Array(items) => {
            for (index, item) in items.iter().enumerate() {
                plain_station_pointers(item, &format!("{at}/{index}"), found);
            }
        }
        _ => {}
    }
}

/// A graph whose run stations write no side and no window is labelled as
/// before 1.2; one side or window on any run station lifts it to 1.2, and
/// a 1.1 reader refuses that payload by the field's name. Every plain
/// station of the every-kind fixture, and every object holding one, is
/// given a window (`0.25` and `-0.0`) and a side in turn, alone: exactly
/// the sections of `SectionsAtStations` and `OpenSectionsAtStations` carry
/// them, and every other station refuses them as unknown fields.
#[test]
fn a_side_or_window_on_any_run_station_and_only_that_selects_one_two() {
    let (graph, _) = fixture();
    let mut payload: Value = serde_json::from_str(&graph.to_json().unwrap()).unwrap();
    assert_eq!(payload["version"], "1.3", "the fixture uses 1.3");
    without_run_seams(&mut payload);
    payload["version"] = json!("1.1");

    let plain = read_json(&payload).unwrap();
    assert_eq!(
        version_of(&plain.to_json().unwrap()),
        "1.1",
        "the fixture's oriented windows"
    );
    assert_eq!(cbor_version_of(&plain.to_cbor().unwrap()), "1.1");

    let mut pointers = Vec::new();
    plain_station_pointers(&payload["graph"]["nodes"], "", &mut pointers);
    let mut carriers = BTreeSet::new();
    let mut others = BTreeSet::new();
    for pointer in &pointers {
        let index: usize = pointer.split('/').nth(1).unwrap().parse().unwrap();
        let kind = kind_of(&payload["graph"]["nodes"][index]);
        for (field, value) in [
            ("seam_window", json!(0.25)),
            ("seam_window", json!(-0.0)),
            ("seam", json!("Incoming")),
        ] {
            let mut lifted = payload.clone();
            lifted["version"] = json!("1.2");
            lifted["graph"]["nodes"]
                .pointer_mut(pointer)
                .unwrap()
                .as_object_mut()
                .unwrap()
                .insert(field.into(), value.clone());
            match read_json(&lifted) {
                Ok(graph) => {
                    carriers.insert(kind.clone());
                    let what = format!("{field} {value} at {pointer}");
                    assert_eq!(version_of(&graph.to_json().unwrap()), "1.2", "{what}");
                    assert_eq!(
                        cbor_version_of(&graph.to_cbor().unwrap()),
                        "1.2",
                        "CBOR: {what}"
                    );
                    lifted["version"] = json!("1.1");
                    refused(
                        &lifted,
                        &WireError::UnknownField {
                            field: field.into(),
                            version: FormatVersion::new(1, 1),
                        },
                    );
                }
                Err(WireError::UnknownField { field: found, .. }) => {
                    assert_eq!(found, field, "{kind} at {pointer}");
                    others.insert(kind.clone());
                }
                Err(other) => panic!("{field} at {pointer}: {other}"),
            }
        }
    }
    assert_eq!(
        carriers,
        BTreeSet::from([
            "SolidOperation::SectionsAtStations".to_owned(),
            "SurfaceRelation::OpenSectionsAtStations".to_owned(),
        ]),
        "every node kind whose run stations carry a side and a window"
    );
    assert_eq!(
        others,
        BTreeSet::from([
            "CurveRelation::OffsetByStations".to_owned(),
            "CurveStation".to_owned(),
            "InstanceAtStation".to_owned(),
            "OrientedCurveStation".to_owned(),
            "SolidOperation::SectionsAtStations".to_owned(),
            "SolidOperation::StationedSpine".to_owned(),
            "SurfaceRelation::OpenSectionsAtStations".to_owned(),
            "SurfaceRelation::SectionedSurface".to_owned(),
        ]),
        "every node kind holding a plain station (a section's own station \
         refuses the fields too: they belong to the section)"
    );
}

/// A one-relation payload along a line under `version`, the relation's
/// node JSON given.
fn relation_payload(version: &str, relation: &Value) -> Value {
    json!({
        "format": FORMAT_NAME,
        "version": version,
        "graph": {
            "nodes": [
                {"Curve3": {"Line": {"origin": [0.0, 0.0, 0.0], "direction": [1.0, 0.0, 0.0]}}},
                relation,
            ],
            "roots": [1],
        }
    })
}

/// `OffsetsAtStations` is a 1.2 variant: written under 1.2 even when its
/// stations write no side or window (which are then omitted), and refused
/// by name in an older payload, as an older reader would refuse it.
#[test]
fn an_offset_at_run_stations_is_a_one_two_variant() {
    let mut builder = GeometryGraphBuilder::new();
    let line = builder.push(GeometryNode::Curve3(line3())).unwrap();
    let offset = builder
        .push(GeometryNode::CurveRelation(
            CurveRelation::OffsetsAtStations {
                basis: line,
                stations: vec![
                    RunStation::new(Station::at(0.0)),
                    RunStation::from(Station::at(2.0)),
                ],
                frame: StationFrame::Section,
            },
        ))
        .unwrap();
    let graph = builder.finish(vec![offset]).unwrap();
    let text = graph.to_json().unwrap();
    assert_eq!(version_of(&text), "1.2");
    assert!(!text.contains("seam"), "{text}");
    let payload: Value = serde_json::from_str(&text).unwrap();
    let relation = payload["graph"]["nodes"][1].clone();
    assert_eq!(
        relation["CurveRelation"]["OffsetsAtStations"]["stations"][0],
        json!({"station": {
            "distance": 0.0,
            "offsets": {"lateral": 0.0, "vertical": 0.0, "longitudinal": 0.0}
        }})
    );
    for version in ["1.0", "1.1"] {
        refused(
            &relation_payload(version, &relation),
            &WireError::UnknownKind {
                kind: "OffsetsAtStations".into(),
                version: FormatVersion::parse(version).unwrap(),
            },
        );
    }
    let read = read_json(&relation_payload("1.2", &relation)).unwrap();
    assert_eq!(read.to_json().unwrap(), text);
}

/// Sides and windows on run stations round-trip bit-exactly through both
/// encodings; `Outgoing` and a `+0.0` window are not written, `-0.0` is;
/// a negative window is refused by validation on read.
#[test]
fn run_station_sides_and_windows_round_trip_and_defaults_are_not_written() {
    let windows = [0.0, -0.0, 1e-3, 0.1 + 0.2];
    let sides = [
        SeamSide::Outgoing,
        SeamSide::Incoming,
        SeamSide::Outgoing,
        SeamSide::Incoming,
    ];
    let mut builder = GeometryGraphBuilder::new();
    let line = builder.push(GeometryNode::Curve3(line3())).unwrap();
    let profile = builder.push(GeometryNode::Profile(rectangle())).unwrap();
    let stations: Vec<RunStation> = (0..4)
        .map(|k| {
            RunStation::new(Station::at(k as f64))
                .with_seam_side(sides[k])
                .with_seam_window(windows[k])
        })
        .collect();
    let sections: Vec<SectionAtStation> = stations
        .iter()
        .map(|station| {
            SectionAtStation::new(profile, station.station)
                .with_seam_side(station.seam)
                .with_seam_window(station.seam_window)
        })
        .collect();
    let offset = builder
        .push(GeometryNode::CurveRelation(
            CurveRelation::OffsetsAtStations {
                basis: line,
                stations,
                frame: StationFrame::Plan,
            },
        ))
        .unwrap();
    let solid = builder
        .push(GeometryNode::SolidOperation(
            SolidOperation::SectionsAtStations {
                directrix: line,
                sections,
                frame: StationFrame::Section,
            },
        ))
        .unwrap();
    let graph = builder.finish(vec![offset, solid]).unwrap();
    let text = graph.to_json().unwrap();
    // Three windows of four written in each run, two sides of four.
    assert_eq!(text.matches("seam_window").count(), 6, "{text}");
    assert_eq!(text.matches("\"seam\":").count(), 4, "{text}");
    let want: Vec<(SeamSide, u64)> = (0..4).map(|k| (sides[k], windows[k].to_bits())).collect();
    for read in [
        GeometryGraph::from_json(&text).unwrap(),
        GeometryGraph::from_cbor(&graph.to_cbor().unwrap()).unwrap(),
    ] {
        let mut runs = 0;
        for (_, node) in read.iter() {
            let got: Vec<(SeamSide, u64)> = match node {
                GeometryNode::CurveRelation(CurveRelation::OffsetsAtStations {
                    stations, ..
                }) => stations
                    .iter()
                    .map(|s| (s.seam, s.seam_window.to_bits()))
                    .collect(),
                GeometryNode::SolidOperation(SolidOperation::SectionsAtStations {
                    sections,
                    ..
                }) => sections
                    .iter()
                    .map(|s| (s.seam, s.seam_window.to_bits()))
                    .collect(),
                _ => continue,
            };
            assert_eq!(got, want);
            runs += 1;
        }
        assert_eq!(runs, 2);
    }
    // A window-free, side-free run writes the graph it wrote before 1.2.
    let mut builder = GeometryGraphBuilder::new();
    let line = builder.push(GeometryNode::Curve3(line3())).unwrap();
    let profile = builder.push(GeometryNode::Profile(rectangle())).unwrap();
    let solid = builder
        .push(GeometryNode::SolidOperation(
            SolidOperation::SectionsAtStations {
                directrix: line,
                sections: vec![
                    SectionAtStation::new(profile, Station::at(0.0)),
                    SectionAtStation::new(profile, Station::at(1.0))
                        .with_seam_side(SeamSide::Outgoing)
                        .with_seam_window(0.0),
                ],
                frame: StationFrame::Section,
            },
        ))
        .unwrap();
    let plain = builder.finish(vec![solid]).unwrap().to_json().unwrap();
    assert_eq!(version_of(&plain), "1.0");
    assert!(!plain.contains("seam"), "{plain}");
    // A malformed window is refused by validation, in either carrier.
    let payload: Value = serde_json::from_str(&text).unwrap();
    for (node, path) in [
        (2, "/CurveRelation/OffsetsAtStations/stations/1/seam_window"),
        (
            3,
            "/SolidOperation/SectionsAtStations/sections/1/seam_window",
        ),
    ] {
        let mut bad = payload.clone();
        *bad["graph"]["nodes"][node].pointer_mut(path).unwrap() = json!(-1e-3);
        assert!(
            matches!(
                read_json(&bad),
                Err(WireError::InvalidGraph {
                    node: Some(n),
                    error: GraphError::InvalidStation { .. },
                }) if n == node
            ),
            "{path}"
        );
    }
}

// ---------------------------------------------------------------------------
// Format 1.3: run stations at their basis's ends (#303)
// ---------------------------------------------------------------------------

/// The node kinds holding a plain station that refuse a field of a run
/// station, as the 1.2 carrier test lists them.
fn plain_station_kinds() -> BTreeSet<String> {
    [
        "CurveRelation::OffsetByStations",
        "CurveRelation::OffsetsAtStations",
        "CurveStation",
        "InstanceAtStation",
        "OrientedCurveStation",
        "SolidOperation::SectionsAtStations",
        "SolidOperation::StationedSpine",
        "SurfaceRelation::OpenSectionsAtStations",
        "SurfaceRelation::SectionedSurface",
    ]
    .into_iter()
    .map(str::to_owned)
    .collect()
}

/// A graph whose run stations all state their distance is labelled as
/// before 1.3; one run station at its basis's start or end lifts it to
/// 1.3, and a 1.2 reader refuses that payload by the field's name. Every
/// plain station of the every-kind fixture, and every object holding one,
/// is given `basis_end` (`"Start"` and `"End"`, its distance set to `0`)
/// in turn, alone: exactly the run stations of `OffsetsAtStations` carry
/// it (where its place in the run is wrong, validation refuses it), and
/// every other station refuses it as an unknown field.
#[test]
fn a_basis_end_on_any_run_station_and_only_that_selects_one_three() {
    let (graph, _) = fixture();
    let mut payload: Value = serde_json::from_str(&graph.to_json().unwrap()).unwrap();
    assert_eq!(payload["version"], "1.3", "the fixture uses 1.3");
    without_basis_ends(&mut payload);
    payload["version"] = json!("1.2");
    assert!(!payload.to_string().contains("basis_end"));

    let plain = read_json(&payload).unwrap();
    assert_eq!(
        version_of(&plain.to_json().unwrap()),
        "1.2",
        "the fixture's run seams"
    );
    assert_eq!(cbor_version_of(&plain.to_cbor().unwrap()), "1.2");

    let mut pointers = Vec::new();
    plain_station_pointers(&payload["graph"]["nodes"], "", &mut pointers);
    let mut carriers = BTreeSet::new();
    let mut lifted_kinds = BTreeSet::new();
    let mut others = BTreeSet::new();
    for pointer in &pointers {
        let index: usize = pointer.split('/').nth(1).unwrap().parse().unwrap();
        let kind = kind_of(&payload["graph"]["nodes"][index]);
        for end in ["Start", "End"] {
            let mut lifted = payload.clone();
            lifted["version"] = json!("1.3");
            let object = lifted["graph"]["nodes"]
                .pointer_mut(pointer)
                .unwrap()
                .as_object_mut()
                .unwrap();
            object.insert("basis_end".into(), json!(end));
            match object.get_mut("station") {
                Some(station) if station.get("distance").is_some() => {
                    station["distance"] = json!(0.0);
                }
                _ => {
                    object.insert("distance".into(), json!(0.0));
                }
            }
            let what = format!("basis_end {end} at {pointer}");
            match read_json(&lifted) {
                Ok(graph) => {
                    carriers.insert(kind.clone());
                    lifted_kinds.insert(kind.clone());
                    assert_eq!(version_of(&graph.to_json().unwrap()), "1.3", "{what}");
                    assert_eq!(
                        cbor_version_of(&graph.to_cbor().unwrap()),
                        "1.3",
                        "CBOR: {what}"
                    );
                    lifted["version"] = json!("1.2");
                    refused(
                        &lifted,
                        &WireError::UnknownField {
                            field: "basis_end".into(),
                            version: FormatVersion::new(1, 2),
                        },
                    );
                }
                Err(WireError::InvalidGraph {
                    error: GraphError::InvalidStation { .. },
                    ..
                }) => {
                    carriers.insert(kind.clone());
                }
                Err(WireError::UnknownField { field, .. }) => {
                    assert_eq!(field, "basis_end", "{what}");
                    others.insert(kind.clone());
                }
                Err(other) => panic!("{what}: {other}"),
            }
        }
    }
    let offsets = BTreeSet::from(["CurveRelation::OffsetsAtStations".to_owned()]);
    assert_eq!(
        carriers, offsets,
        "every node kind whose run stations stand at a basis end"
    );
    assert_eq!(lifted_kinds, offsets, "and each one reads in place");
    assert_eq!(
        others,
        plain_station_kinds(),
        "every node kind holding a plain station (a run station's own \
         station refuses the field too: it belongs to the run station)"
    );
}

/// Run stations at their basis's ends round-trip through both encodings;
/// a stated station writes no `basis_end`, so a run of them stays 1.2; a
/// basis end it does not know is refused by name, one stating a distance
/// by validation, and every one in a 1.2 payload as an unknown field.
#[test]
fn run_station_basis_ends_round_trip_and_a_stated_station_writes_none() {
    let offsets = StationOffsets::new(1.5, -0.25, 0.0);
    let stations = vec![
        RunStation::at_basis_start(offsets),
        RunStation::new(Station::new(4.0, offsets)).with_seam_side(SeamSide::Incoming),
        RunStation::at_basis_end(offsets).with_seam_window(1e-3),
    ];
    let mut builder = GeometryGraphBuilder::new();
    let line = builder.push(GeometryNode::Curve3(line3())).unwrap();
    let offset = builder
        .push(GeometryNode::CurveRelation(
            CurveRelation::OffsetsAtStations {
                basis: line,
                stations: stations.clone(),
                frame: StationFrame::Section,
            },
        ))
        .unwrap();
    let graph = builder.finish(vec![offset]).unwrap();
    let text = graph.to_json().unwrap();
    assert_eq!(version_of(&text), "1.3");
    assert_eq!(text.matches("\"basis_end\":\"Start\"").count(), 1, "{text}");
    assert_eq!(text.matches("\"basis_end\":\"End\"").count(), 1, "{text}");
    assert_eq!(text.matches("basis_end").count(), 2, "{text}");
    for read in [
        GeometryGraph::from_json(&text).unwrap(),
        GeometryGraph::from_cbor(&graph.to_cbor().unwrap()).unwrap(),
    ] {
        let Some(GeometryNode::CurveRelation(CurveRelation::OffsetsAtStations {
            stations: got,
            ..
        })) = read.get(read.roots()[0])
        else {
            panic!("not an offset at run stations");
        };
        assert_eq!(*got, stations);
        assert_eq!(read.to_json().unwrap(), text);
    }
    let payload: Value = serde_json::from_str(&text).unwrap();
    let relation = payload["graph"]["nodes"][1].clone();
    refused(
        &relation_payload("1.2", &relation),
        &WireError::UnknownField {
            field: "basis_end".into(),
            version: FormatVersion::new(1, 2),
        },
    );
    for version in ["1.0", "1.1"] {
        refused(
            &relation_payload(version, &relation),
            &WireError::UnknownKind {
                kind: "OffsetsAtStations".into(),
                version: FormatVersion::parse(version).unwrap(),
            },
        );
    }
    // An unknown end, and an end that states a distance.
    let pointer = "/CurveRelation/OffsetsAtStations/stations/2";
    let mut unknown = relation.clone();
    unknown.pointer_mut(pointer).unwrap()["basis_end"] = json!("Middle");
    refused(
        &relation_payload("1.3", &unknown),
        &WireError::UnknownKind {
            kind: "Middle".into(),
            version: FormatVersion::new(1, 3),
        },
    );
    let mut stated = relation.clone();
    stated.pointer_mut(pointer).unwrap()["station"]["distance"] = json!(9.0);
    assert!(
        matches!(
            read_json(&relation_payload("1.3", &stated)),
            Err(WireError::InvalidGraph {
                node: Some(1),
                error: GraphError::InvalidStation { .. },
            })
        ),
        "an end stating a distance"
    );
    // A run of stated stations writes no basis end: the 1.2 graph.
    let mut builder = GeometryGraphBuilder::new();
    let line = builder.push(GeometryNode::Curve3(line3())).unwrap();
    let offset = builder
        .push(GeometryNode::CurveRelation(
            CurveRelation::OffsetsAtStations {
                basis: line,
                stations: vec![
                    RunStation::new(Station::at(0.0)),
                    RunStation::new(Station::at(2.0)),
                ],
                frame: StationFrame::Section,
            },
        ))
        .unwrap();
    let plain = builder.finish(vec![offset]).unwrap().to_json().unwrap();
    assert_eq!(version_of(&plain), "1.2");
    assert!(!plain.contains("basis_end"), "{plain}");
}
