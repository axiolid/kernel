//! A curved wall less a niche flush with its inner face (#324).
//!
//! A quarter-ring wall (radii `4.8` to `5` about the origin, `0` to `90`
//! degrees, `3` high) less a niche cut from its inner face (radii up to
//! `4.9`, `40` to `50` degrees, `z` from `1` to `2`), whose inner arc is
//! the wall's own (`4.8`), one unit in the last place inside it
//! (`4.8.next_down()`), or well past it (`4.7`). The general
//! boolean refused the first two ("exact boolean whose kept faces do not
//! sew into a solid"), so the deviation report left them unbounded: it
//! read the niche's inner face, whose loops wind clockwise, with its
//! normal turned round where it lies on the wall's. Both are now
//! certified, the niche's inner cylinder read as the wall's within the
//! tolerance and that reading reported, and the exact result measures to
//! the closed-form volume of the annular sectors.

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{BooleanOperator, Frame2, Interval, Point2, Tolerance, Vec2, Vec3};
use axiolid_curve::{Circle2, Curve2, Line2};
use axiolid_measure::exact_properties;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::deviation::COINCIDENT_OPERANDS;
use axiolid_mesh_compile::{DeviationBound, ReferenceExactCompiler, ReferenceMeshCompiler};
use axiolid_model::{
    GeometryGraph, GeometryGraphBuilder, GeometryNode, Instance, NodeId, SolidOperation,
};
use axiolid_profile::{Contour, ContourProfile, Profile, ProfileSegment};

/// An annular sector between radii `inner` and `outer` from angle `from`
/// to `to`, anticlockwise: outer arc, line in, inner arc back, line out,
/// as `ifc-geometry` lowers an `IfcIndexedPolyCurve` of two arcs and two
/// lines.
fn sector(outer: f64, inner: f64, from: f64, to: f64) -> Profile {
    let circle = |radius| {
        Curve2::Circle(Circle2 {
            frame: Frame2 {
                origin: Point2::new(0.0, 0.0),
                x: Vec2::X,
                y: Vec2::Y,
            },
            radius,
        })
    };
    let at = |radius: f64, angle: f64| Point2::new(radius * angle.cos(), radius * angle.sin());
    let line = |a: Point2, b: Point2| ProfileSegment {
        curve: Curve2::Line(Line2 {
            origin: a,
            direction: b - a,
        }),
        domain: Interval::UNIT,
        same_sense: true,
    };
    Profile::Contour(ContourProfile {
        outer: Contour::new(vec![
            ProfileSegment {
                curve: circle(outer),
                domain: Interval::new(from, to),
                same_sense: true,
            },
            line(at(outer, to), at(inner, to)),
            ProfileSegment {
                curve: circle(inner),
                domain: Interval::new(from, to),
                same_sense: false,
            },
            line(at(inner, from), at(outer, from)),
        ]),
        holes: Vec::new(),
    })
}

fn extruded(builder: &mut GeometryGraphBuilder, profile: Profile, z: f64, depth: f64) -> NodeId {
    let profile = builder.push(GeometryNode::Profile(profile)).unwrap();
    let solid = builder
        .push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth,
        }))
        .unwrap();
    builder
        .push(GeometryNode::Instance(Instance {
            source: solid,
            transform: axiolid_core::Transform3::from_translation(Vec3::new(0.0, 0.0, z)),
        }))
        .unwrap()
}

/// The wall less a niche whose inner arc has radius `niche_inner`.
fn wall_less_niche(niche_inner: f64) -> (GeometryGraph, NodeId) {
    let (from, to) = (40.0_f64.to_radians(), 50.0_f64.to_radians());
    let mut builder = GeometryGraphBuilder::new();
    let wall = extruded(
        &mut builder,
        sector(5.0, 4.8, 0.0, std::f64::consts::FRAC_PI_2),
        0.0,
        3.0,
    );
    let niche = extruded(&mut builder, sector(4.9, niche_inner, from, to), 1.0, 1.0);
    let root = builder
        .push(GeometryNode::SolidOperation(SolidOperation::Boolean {
            left: wall,
            right: niche,
            operator: BooleanOperator::Difference,
        }))
        .unwrap();
    (builder.finish(vec![root]).unwrap(), root)
}

/// The wall's volume less the part of the niche inside it.
fn expected_volume(niche_inner: f64) -> f64 {
    let sector = |outer: f64, inner: f64, sweep: f64, height: f64| {
        0.5 * sweep * (outer * outer - inner * inner) * height
    };
    sector(5.0, 4.8, std::f64::consts::FRAC_PI_2, 3.0)
        - sector(4.9, niche_inner.max(4.8), 10.0_f64.to_radians(), 1.0)
}

/// The certified bound of the cut, and how far the operands were read as
/// one surface (`None` where they were not), after checking the exact
/// result's volume.
fn certified(niche_inner: f64) -> (f64, Option<f64>) {
    let options = ExecutionOptions::new(Tolerance::MILLIMETRE);
    let (graph, root) = wall_less_niche(niche_inner);
    let (exact, _) = ReferenceExactCompiler::new()
        .compile_exact_with_report(&graph, root, &options)
        .unwrap_or_else(|e| panic!("niche from {niche_inner}: {e}"));
    let volume = exact_properties(&exact, options.tolerance())
        .expect("measurable")
        .signed_volume;
    let want = expected_volume(niche_inner);
    assert!(
        (volume - want).abs() <= 1e-9 * want,
        "niche from {niche_inner}: expected {want}, got {volume}"
    );
    let compiler = ReferenceMeshCompiler::new(BoolmeshBoolean);
    let (_, report) = compiler
        .compile_mesh_with_deviation(&graph, root, &options)
        .unwrap_or_else(|e| panic!("niche from {niche_inner}: {e}"));
    let bound = report
        .bound
        .unwrap_or_else(|| panic!("niche from {niche_inner}: {report:?}"));
    assert!(bound <= 1e-3, "niche from {niche_inner}: {report:?}");
    let coincident = report
        .contributions
        .iter()
        .find(|c| c.detail == COINCIDENT_OPERANDS)
        .map(|c| match c.bound {
            DeviationBound::Certified(moved) => moved,
            ref other => panic!("niche from {niche_inner}: {other:?}"),
        });
    (bound, coincident)
}

#[test]
fn a_niche_on_the_inner_face_is_certified() {
    let (_, moved) = certified(4.8);
    let moved = moved.expect("the niche's inner cylinder read as the wall's");
    assert!(moved <= 1e-14, "{moved}");
}

#[test]
fn a_niche_one_ulp_inside_the_inner_face_is_certified() {
    let (_, moved) = certified(4.8_f64.next_down());
    let moved = moved.expect("the niche's inner cylinder read as the wall's");
    assert!(moved > 0.0 && moved <= 1e-14, "{moved}");
}

#[test]
fn a_niche_well_inside_is_unchanged() {
    // Certified at 0.99973 mm before #324, with no face read as another's.
    let (bound, moved) = certified(4.7);
    assert_eq!(moved, None);
    assert!((bound - 0.000_999_73).abs() <= 1e-8, "{bound}");
}
