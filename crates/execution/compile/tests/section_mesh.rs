//! Mesh extrusions of structural sections (#193).
//!
//! Each family is compiled through the reference mesh compiler, which
//! refused every `Profile::Section` before. The mesh must be a closed
//! two-manifold (`volume_properties` audits that before it integrates)
//! whose volume is the exact contour's area times the depth: exactly for a
//! section without radii, and within the arcs' chord budget otherwise.

use axiolid_construct::section_lower::{rectangle_contour, section_contour};
use axiolid_contracts::ExecutionOptions;
use axiolid_core::{Tolerance, Transform3, Vec3};
use axiolid_curve::Curve2;
use axiolid_measure::{mesh_distance, volume_properties};
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::ReferenceMeshCompiler;
use axiolid_mesh_compile_contract::MeshCompiler;
use axiolid_model::{GeometryGraphBuilder, GeometryNode, Instance, SolidOperation};
use axiolid_profile::{Contour, ContourProfile, Profile, RectangleProfile, SectionProfile};

const DEPTH: f64 = 2.0;
const CHORD: f64 = 1e-4;

/// The audit's tolerance, finer than the chord budget: a triangle of three
/// neighbouring points on a chorded fillet is only as thick as the chords'
/// sag, and a coarser audit reads it as degenerate.
fn audit() -> Tolerance {
    Tolerance::new(1e-7, 1e-9).expect("a tolerance")
}

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::MILLIMETRE)
        .with_chord_error(CHORD)
        .expect("a positive chord budget")
}

fn extrude(profile: Profile, place: Option<Transform3>) -> TriMesh {
    let mut b = GeometryGraphBuilder::new();
    let profile = b.push(GeometryNode::Profile(profile)).unwrap();
    let mut root = b
        .push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth: DEPTH,
        }))
        .unwrap();
    if let Some(transform) = place {
        root = b
            .push(GeometryNode::Instance(Instance {
                source: root,
                transform,
            }))
            .unwrap();
    }
    let graph = b.finish(vec![root]).unwrap();
    ReferenceMeshCompiler::new(BoolmeshBoolean::new())
        .compile_mesh(&graph, root, &options())
        .expect("a section extrusion meshes")
}

/// Twice the signed area and the arc length of one exact contour ring.
fn ring(contour: &Contour) -> (f64, f64) {
    let (mut twice, mut arcs) = (0.0, 0.0);
    for segment in &contour.segments {
        let (from, to) = if segment.same_sense {
            (segment.domain.start, segment.domain.end)
        } else {
            (segment.domain.end, segment.domain.start)
        };
        match &segment.curve {
            Curve2::Line(line) => {
                let a = line.origin + line.direction * from;
                let b = line.origin + line.direction * to;
                twice += a.perp_dot(b);
            }
            Curve2::Circle(circle) => {
                let point = |t: f64| {
                    let (s, c) = t.sin_cos();
                    circle.frame.origin
                        + circle.frame.x * (circle.radius * c)
                        + circle.frame.y * (circle.radius * s)
                };
                let sweep = (to - from) * circle.frame.x.perp_dot(circle.frame.y).signum();
                twice += point(from).perp_dot(point(to));
                twice += circle.radius * circle.radius * (sweep - sweep.sin());
                arcs += circle.radius * sweep.abs();
            }
            other => panic!("unexpected segment {other:?}"),
        }
    }
    (twice, arcs)
}

/// The exact contour's area, and the total length of its arcs.
fn area_and_arcs(contour: &ContourProfile) -> (f64, f64) {
    let (outer, mut arcs) = ring(&contour.outer);
    let mut area = outer.abs() / 2.0;
    for hole in &contour.holes {
        let (twice, length) = ring(hole);
        area -= twice.abs() / 2.0;
        arcs += length;
    }
    (area, arcs)
}

fn check(name: &str, profile: Profile, contour: &ContourProfile) {
    let mesh = extrude(profile, None);
    let volume = volume_properties(&mesh, audit())
        .unwrap_or_else(|e| panic!("{name}: not a closed solid: {e}"))
        .signed_volume;
    let (area, arcs) = area_and_arcs(contour);
    let exact = area * DEPTH;
    // Each chord strays at most the budget from its arc, so the area
    // differs by at most the arcs' length times the budget.
    let bound = arcs * CHORD * DEPTH + 1e-12 * exact.max(1.0);
    assert!(
        (volume - exact).abs() <= bound,
        "{name}: volume {volume} vs {exact} (bound {bound})"
    );
    if arcs == 0.0 {
        assert!((volume - exact).abs() <= 1e-12 * exact, "{name}: not exact");
    }
}

fn i_section(fillet: Option<f64>, edge: Option<f64>) -> SectionProfile {
    SectionProfile::I {
        depth: 0.3,
        width: 0.3,
        web_thickness: 0.011,
        flange_thickness: 0.019,
        fillet_radius: fillet,
        flange_edge_radius: edge,
        flange_slope: None,
    }
}

#[test]
fn every_section_family_meshes_as_a_closed_prism() {
    let families = [
        ("I, sharp", i_section(None, None)),
        ("I, filleted", i_section(Some(0.027), None)),
        (
            "I, filleted with toe radii",
            i_section(Some(0.027), Some(0.009)),
        ),
        (
            "asymmetric I",
            SectionProfile::AsymmetricI {
                depth: 0.4,
                web_thickness: 0.011,
                bottom_flange_width: 0.3,
                bottom_flange_thickness: 0.019,
                bottom_fillet_radius: Some(0.021),
                bottom_flange_edge_radius: None,
                bottom_flange_slope: None,
                top_flange_width: 0.2,
                top_flange_thickness: Some(0.015),
                top_fillet_radius: Some(0.015),
                top_flange_edge_radius: None,
                top_flange_slope: None,
            },
        ),
        (
            "L",
            SectionProfile::L {
                depth: 0.1,
                width: Some(0.08),
                thickness: 0.008,
                fillet_radius: Some(0.012),
                edge_radius: Some(0.006),
                leg_slope: None,
            },
        ),
        (
            "T",
            SectionProfile::T {
                depth: 0.12,
                flange_width: 0.12,
                web_thickness: 0.011,
                flange_thickness: 0.011,
                fillet_radius: Some(0.011),
                flange_edge_radius: Some(0.005),
                web_edge_radius: Some(0.005),
                web_slope: None,
                flange_slope: None,
            },
        ),
        (
            "U",
            SectionProfile::U {
                depth: 0.3,
                flange_width: 0.1,
                web_thickness: 0.0075,
                flange_thickness: 0.0125,
                fillet_radius: Some(0.015),
                edge_radius: Some(0.006),
                flange_slope: None,
            },
        ),
        (
            "C",
            SectionProfile::C {
                depth: 0.2,
                width: 0.075,
                wall_thickness: 0.003,
                girth: 0.02,
                internal_fillet_radius: Some(0.004),
            },
        ),
        (
            "Z",
            SectionProfile::Z {
                depth: 0.2,
                flange_width: 0.08,
                web_thickness: 0.008,
                flange_thickness: 0.011,
                fillet_radius: Some(0.011),
                edge_radius: Some(0.005),
            },
        ),
        (
            "trapezium",
            SectionProfile::Trapezium {
                bottom_x: 0.4,
                top_x: 0.2,
                y: 0.3,
                top_offset: 0.05,
            },
        ),
    ];
    for (name, section) in families {
        let contour = section_contour(&section).expect("a section lowers");
        check(name, Profile::Section(section), &contour);
    }
}

#[test]
fn rounded_and_hollow_rectangles_mesh_with_their_radii() {
    // Before, the mesh path dropped the corner radii and meshed a sharp
    // box: plausible, and too big.
    for rectangle in [
        RectangleProfile {
            x: 0.2,
            y: 0.1,
            thickness: None,
            outer_radius: Some(0.02),
            inner_radius: None,
        },
        RectangleProfile {
            x: 0.2,
            y: 0.1,
            thickness: Some(0.008),
            outer_radius: Some(0.016),
            inner_radius: Some(0.008),
        },
        RectangleProfile {
            x: 0.2,
            y: 0.1,
            thickness: Some(0.008),
            outer_radius: None,
            inner_radius: None,
        },
    ] {
        let contour = rectangle_contour(&rectangle).expect("a rectangle lowers");
        check(
            &format!("{rectangle:?}"),
            Profile::Rectangle(rectangle),
            &contour,
        );
    }
}

#[test]
fn an_i_beam_resting_on_a_column_touches_it() {
    // A 3 m column, 0.3 x 0.3, standing on z = 0 (extruded 2 m by the
    // helper, so scaled 1.5 along z), and an HEB-like I beam 4 m long along
    // x, its bottom flange on the column's top face at z = 3.
    let column = extrude(
        Profile::Rectangle(RectangleProfile {
            x: 0.3,
            y: 0.3,
            thickness: None,
            outer_radius: None,
            inner_radius: None,
        }),
        Some(Transform3::from_scale(Vec3::new(1.0, 1.0, 1.5))),
    );
    // Profile x (width) along world y, profile y (depth) up, extrusion
    // along world x; the 0.3 deep section's bottom at z = 3.
    let place = |lift: f64| {
        Transform3::from_cols(Vec3::Y, Vec3::Z, Vec3::X, Vec3::new(-1.0, 0.0, 3.15 + lift))
    };
    let beam = extrude(
        Profile::Section(i_section(Some(0.027), None)),
        Some(place(0.0)),
    );
    volume_properties(&beam, audit()).expect("the beam is a closed solid");
    let touching = mesh_distance(&column, &beam).expect("both meshes are usable");
    let d = touching.distance_squared.sqrt();
    assert!(d.is_finite() && d <= 1e-9, "{touching:?}");
    // Lifted 10 mm, the gap is the lift: flat flange over flat top.
    let beam = extrude(
        Profile::Section(i_section(Some(0.027), None)),
        Some(place(0.01)),
    );
    let gap = mesh_distance(&column, &beam).expect("both meshes are usable");
    assert!(
        (gap.distance_squared.sqrt() - 0.01).abs() <= 1e-9,
        "{gap:?}"
    );
}
