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
    use axiolid::mesh_section::{
        ChainEnd, InPlaneSide, RegionLimits, SectionLimits, SectionRequest, SectionRole,
    };
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

    // An open mesh needs the sheet role: this unit square of two triangles is cut
    // by the vertical plane x = 0.5 along a line that ends on its two
    // boundary edges.
    let strip = TriMesh::new(
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ],
        vec![0, 1, 2, 0, 2, 3],
    );
    let vertical = Frame3 {
        origin: Point3::new(0.5, 0.0, 0.0),
        x: Vec3::Y,
        y: Vec3::Z,
        z: Vec3::X,
    };
    let sheet = SectionRequest::new(
        vertical,
        SectionLimits::new(100, 100, 100, 10),
        RegionLimits::new(10, 100, 100),
    );
    assert!(application
        .section_mesh_request(&strip, &sheet, &options)
        .is_err());
    let section = application
        .section_mesh_request(&strip, &sheet.with_role(SectionRole::Sheet), &options)
        .expect("sheet section");
    assert_eq!(section.chains.len(), 1);
    assert_eq!(
        section.chains[0].ends,
        Some([ChainEnd::MeshBoundary, ChainEnd::MeshBoundary])
    );
    assert_eq!(section.evidence.boundary_edges, 4);

    let report = RegionConformanceSuite::run(&ScalarSection::new());
    assert!(report.claimed && report.is_success(), "{report}");
}

#[test]
fn section_request_reports_the_cut_of_a_solid_as_a_closed_region() {
    use axiolid::core::{Frame3, Point2};
    use axiolid::mesh_section::{
        CutRegion, CutRegionLimits, RegionLimits, SectionLimits, SectionRequest,
    };

    let application = ApplicationBuilder::new()
        .with_portable_section()
        .expect("scalar section passes registration conformance")
        .build();
    let options = ExecutionOptions::new(Tolerance::new(1.0e-9, 1.0e-12).expect("valid tolerance"));
    // The plane z = 0.5 cuts the unit cube across its middle.
    let frame = Frame3 {
        origin: Point3::new(0.0, 0.0, 0.5),
        x: Vec3::X,
        y: Vec3::Y,
        z: Vec3::Z,
    };
    let request = SectionRequest::new(
        frame,
        SectionLimits::new(100, 100, 100, 10),
        RegionLimits::new(10, 100, 100),
    );
    let mesh = cube([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let plain = application
        .section_mesh_request(&mesh, &request, &options)
        .expect("plain request");
    assert_eq!(plain.cut_regions, None);

    let cut = application
        .section_mesh_request(
            &mesh,
            &request.with_cut_regions(CutRegionLimits::new(10, 100, 100)),
            &options,
        )
        .expect("cut request");
    assert_eq!(
        cut.cut_regions,
        Some(vec![CutRegion::new(
            vec![
                Point2::new(0.0, 0.0),
                Point2::new(1.0, 0.0),
                Point2::new(1.0, 1.0),
                Point2::new(0.0, 1.0),
            ],
            vec![],
        )])
    );
}

mod drawing {
    use axiolid::application::{Application, ApplicationErrorSource};
    use axiolid::contracts::{GeomError, Operation};
    use axiolid::core::{Frame3, Point2};
    use axiolid::mesh_hidden_lines::{
        BehindSurface, EdgeKind, EdgeSelection, FaceGroup, FaceSource, HiddenLineLimits,
        HiddenLineMesh, HiddenLineRequest, MeshId, MeshRole, Occlusion, RefusedMesh, View,
        Visibility, VisibleFaceLimits, VisibleFaceRequest,
    };

    use super::*;

    fn options() -> ExecutionOptions {
        ExecutionOptions::new(Tolerance::new(1.0e-9, 1.0e-12).expect("valid tolerance"))
    }

    fn plan() -> View {
        View::new(Frame3 {
            origin: Point3::new(0.0, 0.0, 0.0),
            x: Vec3::X,
            y: Vec3::Y,
            z: Vec3::Z,
        })
    }

    fn lines_limits() -> HiddenLineLimits {
        HiddenLineLimits::new(8, 1000, 1000, 100_000, 1000)
    }

    fn faces_limits() -> VisibleFaceLimits {
        VisibleFaceLimits::new(8, 1000, 1000, 100_000, 100, 1000)
    }

    fn square(points: [(f64, f64); 4]) -> Vec<Point2> {
        points.iter().map(|&(x, y)| Point2::new(x, y)).collect()
    }

    /// A horizontal rectangle at height `z`.
    fn rectangle(x0: f64, y0: f64, x1: f64, y1: f64, z: f64) -> TriMesh {
        TriMesh::new(
            vec![
                Point3::new(x0, y0, z),
                Point3::new(x1, y0, z),
                Point3::new(x1, y1, z),
                Point3::new(x0, y1, z),
            ],
            vec![0, 1, 2, 0, 2, 3],
        )
    }

    #[test]
    fn a_pane_over_a_floor_is_seen_through_by_the_lines_beneath_it() {
        // The floor [0,4]^2 at z = 0 and the see-through pane [1,3] x [1,5]
        // at z = 1: the floor's edge y = 4 is seen through the pane from
        // x = 3 to x = 1 and plain elsewhere (ADR 0111).
        let (floor, pane) = (
            rectangle(0.0, 0.0, 4.0, 4.0, 0.0),
            rectangle(1.0, 1.0, 3.0, 5.0, 1.0),
        );
        let glass = [Occlusion::SeeThrough; 2];
        let meshes = [
            HiddenLineMesh::new(MeshId(1), MeshRole::Sheet, &floor),
            HiddenLineMesh::new(MeshId(2), MeshRole::Sheet, &pane).with_occlusion(&glass),
        ];
        let request =
            HiddenLineRequest::new(&meshes, plan(), EdgeSelection::new(None), lines_limits());
        let outcome = Application::portable()
            .expect("portable application")
            .hidden_lines(&request, &options())
            .expect("the pane draws");
        let pieces: Vec<_> = outcome
            .segments
            .iter()
            .filter(|s| s.source.mesh == MeshId(1) && s.source.vertices == [2, 3])
            .map(|s| (s.visibility, s.seen_through, s.start.x, s.end.x))
            .collect();
        let v = Visibility::Visible;
        assert_eq!(
            pieces,
            [(v, 0, 4.0, 3.0), (v, 1, 3.0, 1.0), (v, 0, 1.0, 0.0)]
        );
        assert_eq!(outcome.evidence.see_through_triangles, 2);
        assert_eq!(outcome.evidence.seen_through_segments, 1);
    }

    #[test]
    fn a_pane_over_a_box_lists_the_box_behind_it() {
        let (block, pane) = (
            cube([0.0, 0.0, 0.0], [4.0, 4.0, 1.0]),
            rectangle(1.0, 1.0, 3.0, 5.0, 2.0),
        );
        let glass = [Occlusion::SeeThrough; 2];
        let meshes = [
            HiddenLineMesh::new(MeshId(1), MeshRole::Solid, &block),
            HiddenLineMesh::new(MeshId(2), MeshRole::Sheet, &pane).with_occlusion(&glass),
        ];
        let application = Application::portable().expect("portable application");
        // Nothing may be listed until the caller allows it.
        let refused = VisibleFaceRequest::new(&meshes, plan(), faces_limits());
        let error = application
            .visible_faces(&refused, &options())
            .expect_err("no behind surface allowed");
        let ApplicationErrorSource::HiddenLines(source) = *error.source else {
            panic!("a drawing refusal");
        };
        assert_eq!(
            source.error,
            GeomError::BudgetExceeded {
                resource: "visible-face behind surfaces"
            }
        );
        let request =
            VisibleFaceRequest::new(&meshes, plan(), faces_limits().with_max_behind_surfaces(4));
        let outcome = application
            .visible_faces(&request, &options())
            .expect("the pane resolves");
        let stacks: Vec<_> = outcome
            .regions
            .iter()
            .map(|r| (r.mesh, r.behind.clone()))
            .collect();
        let box_top = BehindSurface::new(MeshId(1), FaceSource::Group(FaceGroup(0)), 1);
        assert_eq!(
            stacks,
            [
                (MeshId(1), vec![]),
                (MeshId(2), vec![box_top]),
                (MeshId(2), vec![])
            ]
        );
        assert_eq!(outcome.evidence.behind_surfaces, 1);
    }

    fn triangle() -> TriMesh {
        TriMesh::new(
            vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
            ],
            vec![0, 1, 2],
        )
    }

    #[test]
    fn a_perspective_cube_draws_its_front_ring_visible_and_the_rest_hidden() {
        // The box [-1,1]^2 x [-4,-2] seen from the eye at the origin: the
        // front ring projects to +-1/2 and the back ring to +-1/4.
        let mesh = cube([-1.0, -1.0, -4.0], [1.0, 1.0, -2.0]);
        let meshes = [HiddenLineMesh::new(MeshId(7), MeshRole::Solid, &mesh)];
        let view = View::perspective(
            Frame3 {
                origin: Point3::new(0.0, 0.0, 0.0),
                x: Vec3::X,
                y: Vec3::Y,
                z: Vec3::Z,
            },
            -1.0,
        )
        .expect("a negative near");
        let request =
            HiddenLineRequest::new(&meshes, view, EdgeSelection::new(Some(0.5)), lines_limits());
        let application = Application::portable().expect("portable application");
        let outcome = application
            .hidden_lines(&request, &options())
            .expect("the perspective cube draws");
        assert_eq!(outcome.view, view);
        let count = |visibility| {
            outcome
                .segments
                .iter()
                .filter(|s| s.visibility == visibility)
                .count()
        };
        assert_eq!(
            (count(Visibility::Visible), count(Visibility::Hidden)),
            (4, 8)
        );
        for segment in outcome
            .segments
            .iter()
            .filter(|s| s.visibility == Visibility::Visible)
        {
            for p in [segment.start, segment.end] {
                assert_eq!((p.x.abs(), p.y.abs()), (0.5, 0.5));
            }
        }

        // Visible faces still refuse perspective, by routing.
        let faces = VisibleFaceRequest::new(&meshes, view, faces_limits());
        let error = application
            .visible_faces(&faces, &options())
            .expect_err("no visible-face provider draws perspective");
        let ApplicationErrorSource::HiddenLines(source) = *error.source else {
            panic!("a drawing refusal");
        };
        assert!(
            matches!(source.error, GeomError::Unsupported { .. }),
            "{source}"
        );
    }

    #[test]
    fn the_plan_cube_draws_four_hidden_and_four_visible_silhouettes() {
        let mesh = cube([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        let meshes = [HiddenLineMesh::new(MeshId(7), MeshRole::Solid, &mesh)];
        let request = HiddenLineRequest::new(
            &meshes,
            plan(),
            EdgeSelection::new(Some(0.5)),
            lines_limits(),
        );
        let outcome = Application::portable()
            .expect("portable application")
            .hidden_lines(&request, &options())
            .expect("the cube draws");

        let (h, v) = (Visibility::Hidden, Visibility::Visible);
        let expected = [
            ([0, 1], h, (0.0, 0.0), (2.0, 0.0)),
            ([0, 3], h, (0.0, 0.0), (0.0, 2.0)),
            ([1, 2], h, (2.0, 0.0), (2.0, 2.0)),
            ([2, 3], h, (2.0, 2.0), (0.0, 2.0)),
            ([4, 5], v, (0.0, 0.0), (2.0, 0.0)),
            ([4, 7], v, (0.0, 0.0), (0.0, 2.0)),
            ([5, 6], v, (2.0, 0.0), (2.0, 2.0)),
            ([6, 7], v, (2.0, 2.0), (0.0, 2.0)),
        ];
        assert_eq!(outcome.segments.len(), expected.len());
        for (segment, (vertices, visibility, start, end)) in outcome.segments.iter().zip(expected) {
            assert_eq!(segment.source.mesh, MeshId(7));
            assert_eq!(segment.source.vertices, vertices);
            assert_eq!(segment.kind, EdgeKind::Silhouette);
            assert_eq!(segment.visibility, visibility);
            assert_eq!(segment.start, Point2::new(start.0, start.1));
            assert_eq!(segment.end, Point2::new(end.0, end.1));
        }
        let evidence = outcome.evidence;
        assert_eq!(evidence.selected_edges, 12);
        assert_eq!(evidence.point_edges, 4);
        assert_eq!(evidence.clipped_edges, 0);
        assert_eq!(evidence.collapsed_segments, 0);
        assert_eq!(evidence.visible_segments, 4);
        assert_eq!(evidence.hidden_segments, 4);
    }

    #[test]
    fn the_plan_cube_cut_by_a_near_plane_has_one_cap_region() {
        let mesh = cube([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        let meshes = [HiddenLineMesh::new(MeshId(7), MeshRole::Solid, &mesh)];
        let request = VisibleFaceRequest::new(&meshes, plan().with_near(1.0), faces_limits());
        let outcome = Application::portable()
            .expect("portable application")
            .visible_faces(&request, &options())
            .expect("the cube cut draws");

        assert_eq!(outcome.regions.len(), 1);
        let region = &outcome.regions[0];
        assert_eq!(region.source, FaceSource::Cap);
        assert_eq!(
            region.outer,
            square([(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)])
        );
        assert!(region.holes.is_empty());
        assert!(region.triangles.is_empty());
    }

    #[test]
    fn the_plan_cube_with_groups_has_one_region_of_its_top_group() {
        let mesh = cube([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        let groups: Vec<FaceGroup> = [0, 0, 2, 2, 4, 4, 6, 6, 8, 8, 10, 10]
            .into_iter()
            .map(FaceGroup)
            .collect();
        let meshes = [HiddenLineMesh::new(MeshId(7), MeshRole::Solid, &mesh).with_groups(&groups)];
        let request = VisibleFaceRequest::new(&meshes, plan(), faces_limits());
        let outcome = Application::portable()
            .expect("portable application")
            .visible_faces(&request, &options())
            .expect("the cube draws");

        assert_eq!(outcome.regions.len(), 1);
        let region = &outcome.regions[0];
        assert_eq!(region.source, FaceSource::Group(FaceGroup(2)));
        assert_eq!(
            region.outer,
            square([(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)])
        );
        assert!(region.holes.is_empty());
        assert_eq!(region.triangles, [2, 3]);
    }

    #[test]
    fn a_refused_mesh_keeps_its_id_in_the_error_source() {
        let mesh = triangle();
        let meshes = [HiddenLineMesh::new(MeshId(9), MeshRole::Solid, &mesh)];
        let request = HiddenLineRequest::new(
            &meshes,
            plan(),
            EdgeSelection::new(Some(0.5)),
            lines_limits(),
        );
        let error = Application::portable()
            .expect("portable application")
            .hidden_lines(&request, &options())
            .expect_err("a lone triangle is not a solid");

        assert_eq!(error.context.operation, Operation::MeshHiddenLines);
        assert_eq!(error.context.provider.as_str(), "scalar-hidden-lines");
        let ApplicationErrorSource::HiddenLines(source) = error.source.as_ref() else {
            panic!("expected a hidden-line refusal, got {:?}", error.source);
        };
        assert_eq!(source.mesh, Some(MeshId(9)));
        assert!(matches!(source.error, GeomError::NotManifold(_)));
    }

    #[test]
    fn a_refused_face_request_keeps_its_operation_provider_and_mesh() {
        let mesh = triangle();
        let meshes = [HiddenLineMesh::new(MeshId(9), MeshRole::Solid, &mesh)];
        let request = VisibleFaceRequest::new(&meshes, plan(), faces_limits());
        let error = Application::portable()
            .expect("portable application")
            .visible_faces(&request, &options())
            .expect_err("a lone triangle is not a solid");

        assert_eq!(error.context.operation, Operation::MeshVisibleFaces);
        assert_eq!(error.context.provider.as_str(), "scalar-visible-faces");
        let ApplicationErrorSource::HiddenLines(source) = error.source.as_ref() else {
            panic!("expected a refusal, got {:?}", error.source);
        };
        assert_eq!(source.mesh, Some(MeshId(9)));
        assert!(matches!(source.error, GeomError::NotManifold(_)));
    }

    #[test]
    fn a_refused_audit_names_its_operation_and_provider() {
        let mesh = triangle();
        let meshes = [HiddenLineMesh::new(MeshId(9), MeshRole::Solid, &mesh)];
        let request = HiddenLineRequest::new(
            &meshes,
            plan(),
            EdgeSelection::new(Some(0.5)),
            HiddenLineLimits::new(0, 1000, 1000, 100_000, 1000),
        );
        let error = Application::portable()
            .expect("portable application")
            .audit_hidden_lines(&request, &options())
            .expect_err("the request exceeds its mesh limit");

        assert_eq!(error.context.operation, Operation::MeshHiddenLines);
        assert_eq!(error.context.provider.as_str(), "scalar-hidden-lines");
    }

    #[test]
    fn the_audit_lists_every_refused_mesh_in_request_order() {
        let solid = cube([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        let lone = triangle();
        let meshes = [
            HiddenLineMesh::new(MeshId(1), MeshRole::Solid, &solid),
            HiddenLineMesh::new(MeshId(9), MeshRole::Solid, &lone),
        ];
        let request = HiddenLineRequest::new(
            &meshes,
            plan(),
            EdgeSelection::new(Some(0.5)),
            lines_limits(),
        );
        let refused = Application::portable()
            .expect("portable application")
            .audit_hidden_lines(&request, &options())
            .expect("the request itself is valid");

        assert_eq!(refused.len(), 1);
        assert_eq!(
            refused[0],
            RefusedMesh::new(MeshId(9), refused[0].error.clone())
        );
        assert!(matches!(refused[0].error, GeomError::NotManifold(_)));
    }

    #[test]
    fn an_application_without_a_provider_refuses_as_unsupported() {
        let mesh = cube([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        let meshes = [HiddenLineMesh::new(MeshId(7), MeshRole::Solid, &mesh)];
        let request = HiddenLineRequest::new(
            &meshes,
            plan(),
            EdgeSelection::new(Some(0.5)),
            lines_limits(),
        );
        let error = ApplicationBuilder::new()
            .build()
            .hidden_lines(&request, &options())
            .expect_err("no provider is selected");

        assert_eq!(error.context.provider.as_str(), "axiolid-application");
        let ApplicationErrorSource::HiddenLines(source) = error.source.as_ref() else {
            panic!("expected a hidden-line refusal, got {:?}", error.source);
        };
        assert!(matches!(
            source.error,
            GeomError::Unsupported { backend, .. } if backend.as_str() == "mesh-hidden-lines-registry"
        ));
    }

    #[test]
    fn the_descriptor_advertises_both_drawing_providers() {
        let application = Application::portable().expect("portable application");
        let descriptor = application.descriptor();
        for (id, provider) in [
            (capability_ids::MESH_HIDDEN_LINES, "scalar-hidden-lines"),
            (capability_ids::MESH_VISIBLE_FACES, "scalar-visible-faces"),
        ] {
            let capability = descriptor
                .require(CapabilityRequirement {
                    id,
                    output: Representation::Profile2d,
                    exactness: Exactness::ToleranceBounded,
                    deterministic: true,
                })
                .expect("drawing capability advertised");
            assert_eq!(capability.provider.as_str(), provider);
            assert_eq!(capability.exactness, Exactness::ToleranceBounded);
            assert_eq!(capability.output, Representation::Profile2d);
        }
        for feature in [
            "mesh-hidden-lines",
            "dispatch-mesh-hidden-lines",
            "hidden-lines-provider",
        ] {
            assert!(descriptor.enabled_features.contains(&feature), "{feature}");
        }
        let bare = ApplicationBuilder::new().build();
        assert!(bare
            .descriptor()
            .capabilities
            .iter()
            .all(|capability| capability.id != capability_ids::MESH_HIDDEN_LINES));
    }
}
