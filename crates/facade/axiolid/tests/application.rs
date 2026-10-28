#![cfg(feature = "application")]

use axiolid::application::ApplicationBuilder;
use axiolid::contracts::{
    capability_ids, CapabilityRequirement, Exactness, ExecutionOptions, Representation,
};
use axiolid::core::{BooleanOperator, Point3, Ray3, Tolerance, Vec3};
use axiolid::mesh::TriMesh;
use axiolid::profile::{Profile, RectangleProfile};

fn cube(min: [f64; 3], max: [f64; 3]) -> TriMesh {
    let [x0, y0, z0] = min;
    let [x1, y1, z1] = max;
    let positions = vec![
        Point3::new(x0, y0, z0),
        Point3::new(x1, y0, z0),
        Point3::new(x1, y1, z0),
        Point3::new(x0, y1, z0),
        Point3::new(x0, y0, z1),
        Point3::new(x1, y0, z1),
        Point3::new(x1, y1, z1),
        Point3::new(x0, y1, z1),
    ];
    let indices = vec![
        0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5, 2, 3, 7, 2, 7, 6,
        3, 0, 4, 3, 4, 7,
    ];
    TriMesh::new(positions, indices)
}

#[test]
fn portable_application_executes_the_reference_workflows() {
    let application = ApplicationBuilder::new()
        .with_portable_boolean()
        .expect("boolmesh passes registration conformance")
        .with_portable_section()
        .expect("scalar section passes registration conformance")
        .build();
    let tolerance = Tolerance::new(1.0e-9, 1.0e-12).expect("valid tolerance");
    let options = ExecutionOptions::new(tolerance);
    let subject = cube([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
    let tool = cube([1.0, 1.0, 0.5], [3.0, 3.0, 1.5]);

    let health = application
        .validate_mesh(&subject, tolerance)
        .expect("mesh audit");
    assert!(health.is_closed_two_manifold());

    let measurements = application
        .measure_mesh(&subject, tolerance)
        .expect("measurement");
    assert!((measurements.volume.signed_volume - 8.0).abs() < 1.0e-10);

    let difference = application
        .boolean(&subject, &tool, BooleanOperator::Difference, &options)
        .expect("boolean difference");
    assert!(!difference.mesh.indices.is_empty());

    let batched = application
        .subtract_many(&subject, std::slice::from_ref(&tool), &options)
        .expect("batched subtraction");
    assert!(!batched.mesh.indices.is_empty());

    let ray = Ray3 {
        origin: Point3::new(-1.0, 1.0, 1.0),
        direction: Vec3::X,
    };
    let hit = application
        .nearest_mesh_hit(&subject, &ray, tolerance)
        .expect("spatial query")
        .expect("ray intersects cube");
    assert!((hit.point.x - 0.0).abs() < 1.0e-10);

    let profile = Profile::Rectangle(RectangleProfile {
        x: 2.0,
        y: 3.0,
        thickness: None,
        outer_radius: None,
        inner_radius: None,
    });
    let exact = application
        .extrude_profile_exact(&profile, Vec3::Z, 4.0, tolerance)
        .expect("exact rectangle extrusion");
    assert!(!exact.topology().faces().is_empty());
}

#[test]
fn descriptor_reports_selected_providers_and_exactness() {
    let application = ApplicationBuilder::new()
        .with_portable_boolean()
        .expect("conformant provider")
        .with_portable_section()
        .expect("scalar section passes registration conformance")
        .build();
    let descriptor = application.descriptor();

    let mesh_boolean = descriptor
        .require(CapabilityRequirement {
            id: capability_ids::MESH_BOOLEAN,
            output: Representation::TriangleMesh,
            exactness: Exactness::ToleranceBounded,
            deterministic: true,
        })
        .expect("mesh boolean advertised");
    assert_eq!(mesh_boolean.provider.as_str(), "boolmesh");

    descriptor
        .require(CapabilityRequirement {
            id: capability_ids::EXACT_EXTRUDE,
            output: Representation::ExactBrep,
            exactness: Exactness::Exact,
            deterministic: true,
        })
        .expect("exact extrusion advertised");
}

#[test]
fn section_request_reports_a_face_in_the_plane_as_a_region() {
    use axiolid::core::{Frame3, Point2};
    use axiolid::mesh_section::conformance::RegionConformanceSuite;
    use axiolid::mesh_section::{InPlaneSide, RegionLimits, SectionLimits, SectionRequest};
    use axiolid_reference::ScalarSection;

    let application = ApplicationBuilder::new()
        .with_portable_section()
        .expect("scalar section passes registration conformance")
        .build();
    let options = ExecutionOptions::new(Tolerance::new(1.0e-9, 1.0e-12).expect("valid tolerance"));
    // The plane z = 1 contains the top face of the unit cube.
    let frame = Frame3 {
        origin: Point3::new(0.0, 0.0, 1.0),
        x: Vec3::X,
        y: Vec3::Y,
        z: Vec3::Z,
    };
    let request = SectionRequest::new(
        frame,
        SectionLimits::new(100, 100, 100, 10),
        RegionLimits::new(10, 100, 100),
    );
    let section = application
        .section_mesh_request(&cube([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]), &request, &options)
        .expect("section request");
    assert!(section.chains.is_empty());
    assert_eq!(section.regions.len(), 1);
    assert_eq!(section.regions[0].side, InPlaneSide::Below);
    assert_eq!(
        section.regions[0].outer,
        vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(0.0, 1.0),
        ]
    );

    let report = RegionConformanceSuite::run(&ScalarSection::new());
    assert!(report.claimed && report.is_success(), "{report}");
}
