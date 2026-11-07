use std::ptr::null;

use axiolid_capi::*;

fn context_with(config: AxiolidContextConfig) -> AxiolidContextHandle {
    unsafe {
        let mut handle = AxiolidContextHandle::INVALID;
        assert_eq!(
            axiolid_v0_4_context_create(&config, &mut handle),
            AxiolidStatus::Ok
        );
        handle
    }
}

fn context() -> AxiolidContextHandle {
    context_with(AxiolidContextConfig::default())
}

fn import(context: AxiolidContextHandle, positions: &[f64], indices: &[u32]) -> AxiolidMeshHandle {
    unsafe {
        let mut mesh = AxiolidMeshHandle::INVALID;
        assert_eq!(
            axiolid_v0_4_mesh_import(
                context,
                positions.as_ptr(),
                positions.len() / 3,
                indices.as_ptr(),
                indices.len() / 3,
                &mut mesh,
            ),
            AxiolidStatus::Ok
        );
        mesh
    }
}

fn cube(context: AxiolidContextHandle) -> AxiolidMeshHandle {
    import(
        context,
        &[
            0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 2.0, 2.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0, 2.0, 0.0,
            2.0, 2.0, 2.0, 2.0, 0.0, 2.0, 2.0,
        ],
        &[
            0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5, 2, 3, 7, 2, 7,
            6, 3, 0, 4, 3, 4, 7,
        ],
    )
}

fn triangle(context: AxiolidContextHandle) -> AxiolidMeshHandle {
    import(
        context,
        &[0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0],
        &[0, 1, 2],
    )
}

/// The unit square as two triangles sharing the diagonal 0-2.
fn square(context: AxiolidContextHandle) -> AxiolidMeshHandle {
    import(
        context,
        &[0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0],
        &[0, 1, 2, 0, 2, 3],
    )
}

fn entry(mesh: AxiolidMeshHandle, id: u64, role: i32) -> AxiolidViewMesh {
    AxiolidViewMesh {
        mesh,
        id,
        groups: null(),
        group_count: 0,
        role,
        group_mode: AXIOLID_GROUP_MODE_SINGLE,
    }
}

fn plan() -> AxiolidView {
    AxiolidView {
        frame: AxiolidViewFrame {
            origin: [0.0; 3],
            x: [1.0, 0.0, 0.0],
            y: [0.0, 1.0, 0.0],
            z: [0.0, 0.0, 1.0],
        },
        near: 0.0,
        far: 0.0,
        projection: AXIOLID_PROJECTION_ORTHOGRAPHIC,
        has_near: 0,
        has_far: 0,
        reserved: [0; 2],
    }
}

fn all_edges() -> AxiolidEdgeSelection {
    AxiolidEdgeSelection {
        crease_cosine: 0.5,
        has_crease: 1,
        boundary: 1,
        silhouette: 1,
        hidden: 1,
        group_boundary: 0,
        reserved: [0; 3],
    }
}

fn limits() -> AxiolidHiddenLineLimits {
    AxiolidHiddenLineLimits {
        max_meshes: 8,
        max_source_vertices: 1000,
        max_source_triangles: 1000,
        max_candidate_pairs: 100_000,
        max_output_segments: 1000,
    }
}

const TOLERANCE: AxiolidTolerance = AxiolidTolerance {
    linear: 1.0e-9,
    angular: 1.0e-12,
};

type Change = Box<dyn Fn(&mut Call)>;

struct Call {
    meshes: Vec<AxiolidViewMesh>,
    view: AxiolidView,
    edges: AxiolidEdgeSelection,
    region: Option<AxiolidImageBox>,
    limits: AxiolidHiddenLineLimits,
    tolerance: AxiolidTolerance,
}

impl Call {
    fn new(meshes: Vec<AxiolidViewMesh>) -> Self {
        Self {
            meshes,
            view: plan(),
            edges: all_edges(),
            region: None,
            limits: limits(),
            tolerance: TOLERANCE,
        }
    }

    fn run(&self, context: AxiolidContextHandle) -> (AxiolidStatus, AxiolidResultHandle) {
        let mut result = AxiolidResultHandle::INVALID;
        let region = self.region.as_ref().map_or(null(), |region| region);
        let status = unsafe {
            axiolid_v0_4_hidden_lines(
                context,
                self.meshes.as_ptr(),
                self.meshes.len(),
                &self.view,
                &self.edges,
                region,
                &self.limits,
                self.tolerance,
                &mut result,
            )
        };
        (status, result)
    }
}

fn segments(
    context: AxiolidContextHandle,
    result: AxiolidResultHandle,
) -> Vec<AxiolidHiddenLineSegment> {
    unsafe {
        let mut count = 0_usize;
        assert_eq!(
            axiolid_v0_4_hidden_lines_count(context, result, &mut count),
            AxiolidStatus::Ok
        );
        let mut out = vec![
            AxiolidHiddenLineSegment {
                mesh_id: 0,
                vertices: [0; 2],
                start: [0.0; 2],
                end: [0.0; 2],
                kind: AxiolidEdgeKind::Unknown,
                visibility: AxiolidVisibility::Unknown,
            };
            count
        ];
        assert_eq!(
            axiolid_v0_4_hidden_lines_copy(context, result, out.as_mut_ptr(), count),
            AxiolidStatus::Ok
        );
        out
    }
}

fn last_error_mesh(context: AxiolidContextHandle) -> (AxiolidStatus, u8, u64) {
    let (mut has, mut id) = (9_u8, 9_u64);
    let status = unsafe { axiolid_v0_4_context_last_error_mesh(context, &mut has, &mut id) };
    (status, has, id)
}

#[test]
fn the_plan_cube_draws_four_hidden_and_four_visible_silhouettes() {
    unsafe {
        let context = context();
        let cube = cube(context);
        let (status, result) =
            Call::new(vec![entry(cube, 7, AXIOLID_MESH_ROLE_SOLID)]).run(context);
        assert_eq!(status, AxiolidStatus::Ok);

        let mut kind = AxiolidGeometryKind::Unknown;
        assert_eq!(
            axiolid_v0_4_result_kind(context, result, &mut kind),
            AxiolidStatus::Ok
        );
        assert_eq!(kind, AxiolidGeometryKind::HiddenLines);
        let mut taken = AxiolidMeshHandle::INVALID;
        assert_eq!(
            axiolid_v0_4_result_take_mesh(context, result, &mut taken),
            AxiolidStatus::WrongResultKind
        );

        let (h, v) = (AxiolidVisibility::Hidden, AxiolidVisibility::Visible);
        let expected = [
            ([0, 1], h, [0.0, 0.0], [2.0, 0.0]),
            ([0, 3], h, [0.0, 0.0], [0.0, 2.0]),
            ([1, 2], h, [2.0, 0.0], [2.0, 2.0]),
            ([2, 3], h, [2.0, 2.0], [0.0, 2.0]),
            ([4, 5], v, [0.0, 0.0], [2.0, 0.0]),
            ([4, 7], v, [0.0, 0.0], [0.0, 2.0]),
            ([5, 6], v, [2.0, 0.0], [2.0, 2.0]),
            ([6, 7], v, [2.0, 2.0], [0.0, 2.0]),
        ];
        let drawn = segments(context, result);
        assert_eq!(drawn.len(), expected.len());
        for (segment, (vertices, visibility, start, end)) in drawn.iter().zip(expected) {
            assert_eq!(segment.mesh_id, 7);
            assert_eq!(segment.vertices, vertices);
            assert_eq!(segment.kind, AxiolidEdgeKind::Silhouette);
            assert_eq!(segment.visibility, visibility);
            assert_eq!((segment.start, segment.end), (start, end));
        }

        let mut evidence = std::mem::zeroed::<AxiolidHiddenLineEvidence>();
        assert_eq!(
            axiolid_v0_4_hidden_lines_evidence(context, result, &mut evidence),
            AxiolidStatus::Ok
        );
        assert_eq!(
            (
                evidence.source_meshes,
                evidence.source_triangles,
                evidence.selected_edges,
                evidence.point_edges,
                evidence.clipped_edges,
                evidence.collapsed_segments,
                evidence.visible_segments,
                evidence.hidden_segments,
            ),
            (1, 12, 12, 4, 0, 0, 4, 4)
        );

        // A short buffer is refused and nothing is written.
        let sentinel = AxiolidHiddenLineSegment {
            mesh_id: 99,
            vertices: [0; 2],
            start: [0.0; 2],
            end: [0.0; 2],
            kind: AxiolidEdgeKind::Unknown,
            visibility: AxiolidVisibility::Unknown,
        };
        let mut short = vec![sentinel; 7];
        assert_eq!(
            axiolid_v0_4_hidden_lines_copy(context, result, short.as_mut_ptr(), 7),
            AxiolidStatus::BufferTooSmall
        );
        assert!(short.iter().all(|value| *value == sentinel));

        assert_eq!(
            axiolid_v0_4_result_destroy(context, result),
            AxiolidStatus::Ok
        );
        assert_eq!(axiolid_v0_4_mesh_destroy(context, cube), AxiolidStatus::Ok);
        let (mut meshes, mut results) = (1_usize, 1_usize);
        assert_eq!(
            axiolid_v0_4_context_live_object_counts(context, &mut meshes, &mut results),
            AxiolidStatus::Ok
        );
        assert_eq!((meshes, results), (0, 0));
        assert_eq!(axiolid_v0_4_context_destroy(context), AxiolidStatus::Ok);
    }
}

#[test]
fn a_depth_bound_is_read_only_when_its_flag_is_set() {
    let context = context();
    let cube = cube(context);
    let mut call = Call::new(vec![entry(cube, 1, AXIOLID_MESH_ROLE_SOLID)]);
    // The value is ignored while the flag is clear.
    call.view.near = -5.0;
    call.view.far = 5.0;
    let (status, result) = call.run(context);
    assert_eq!(status, AxiolidStatus::Ok);
    assert_eq!(segments(context, result).len(), 8);

    // Without hidden pieces only the four visible silhouettes remain.
    let mut visible_only = Call::new(vec![entry(cube, 1, AXIOLID_MESH_ROLE_SOLID)]);
    visible_only.edges.hidden = 0;
    let (status, result) = visible_only.run(context);
    assert_eq!(status, AxiolidStatus::Ok);
    let drawn = segments(context, result);
    assert_eq!(drawn.len(), 4);
    assert!(drawn
        .iter()
        .all(|segment| segment.visibility == AxiolidVisibility::Visible));

    // Everything lies at W >= 0, so a near plane at -5 keeps nothing.
    call.view.has_near = 1;
    let (status, result) = call.run(context);
    assert_eq!(status, AxiolidStatus::Ok);
    assert!(segments(context, result).is_empty());
}

#[test]
fn a_region_crops_the_image_and_null_is_no_crop() {
    let context = context();
    let cube = cube(context);
    let mut call = Call::new(vec![entry(cube, 1, AXIOLID_MESH_ROLE_SOLID)]);
    call.region = Some(AxiolidImageBox {
        min: [10.0, 10.0],
        max: [11.0, 11.0],
    });
    let (status, result) = call.run(context);
    assert_eq!(status, AxiolidStatus::Ok);
    assert!(segments(context, result).is_empty());
}

#[test]
fn per_triangle_groups_draw_the_group_boundary() {
    let context = context();
    let square = square(context);
    let groups = [0_u32, 1];
    let mut mesh = entry(square, 4, AXIOLID_MESH_ROLE_SHEET);
    mesh.groups = groups.as_ptr();
    mesh.group_count = groups.len();
    mesh.group_mode = AXIOLID_GROUP_MODE_PER_TRIANGLE;
    let mut call = Call::new(vec![mesh]);
    call.edges.group_boundary = 1;
    let (status, result) = call.run(context);
    assert_eq!(status, AxiolidStatus::Ok);
    let drawn = segments(context, result);
    assert_eq!(drawn.len(), 5);
    let diagonal = drawn
        .iter()
        .find(|segment| segment.vertices == [0, 2])
        .expect("the diagonal is drawn");
    assert_eq!(diagonal.kind, AxiolidEdgeKind::GroupBoundary);
    assert!(drawn
        .iter()
        .filter(|segment| segment.vertices != [0, 2])
        .all(|segment| segment.kind == AxiolidEdgeKind::Boundary));
    let mut evidence = unsafe { std::mem::zeroed::<AxiolidHiddenLineEvidence>() };
    assert_eq!(
        unsafe { axiolid_v0_4_hidden_lines_evidence(context, result, &mut evidence) },
        AxiolidStatus::Ok
    );
    assert_eq!(
        (evidence.visible_segments, evidence.hidden_segments),
        (5, 0)
    );
    // Every escalation settles one examined pair, so there are more pairs.
    assert!(evidence.exact_escalations > 0);
    assert!(evidence.candidate_pairs > evidence.exact_escalations);

    // Deselected kinds are not drawn: without boundary edges only the diagonal remains.
    call.edges.boundary = 0;
    let (status, result) = call.run(context);
    assert_eq!(status, AxiolidStatus::Ok);
    let drawn = segments(context, result);
    assert_eq!(drawn.len(), 1);
    assert_eq!(drawn[0].vertices, [0, 2]);
    call.edges.boundary = 1;

    // Planar groups derive one patch: the diagonal is no group boundary.
    call.meshes[0] = entry(square, 4, AXIOLID_MESH_ROLE_SHEET);
    call.meshes[0].group_mode = AXIOLID_GROUP_MODE_PLANAR;
    let (status, result) = call.run(context);
    assert_eq!(status, AxiolidStatus::Ok);
    assert_eq!(segments(context, result).len(), 4);

    // Two triangles folded along an edge are two planar groups; with a single
    // group the fold is no group boundary.
    let fold = import(
        context,
        &[0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0],
        &[0, 1, 2, 1, 3, 2],
    );
    let mut folded = Call::new(vec![entry(fold, 5, AXIOLID_MESH_ROLE_SHEET)]);
    folded.edges.has_crease = 0;
    folded.edges.group_boundary = 1;
    let (status, result) = folded.run(context);
    assert_eq!(status, AxiolidStatus::Ok);
    assert_eq!(segments(context, result).len(), 4);
    folded.meshes[0].group_mode = AXIOLID_GROUP_MODE_PLANAR;
    let (status, result) = folded.run(context);
    assert_eq!(status, AxiolidStatus::Ok);
    let drawn = segments(context, result);
    assert_eq!(drawn.len(), 5);
    let shared = drawn.iter().find(|segment| segment.vertices == [1, 2]);
    assert_eq!(
        shared.map(|segment| segment.kind),
        Some(AxiolidEdgeKind::GroupBoundary)
    );

    // A fold is a crease once the threshold admits its angle.
    let mut creased = Call::new(vec![entry(fold, 5, AXIOLID_MESH_ROLE_SHEET)]);
    creased.edges.crease_cosine = 0.9;
    let (status, result) = creased.run(context);
    assert_eq!(status, AxiolidStatus::Ok);
    let drawn = segments(context, result);
    let shared = drawn.iter().find(|segment| segment.vertices == [1, 2]);
    assert_eq!(
        shared.map(|segment| segment.kind),
        Some(AxiolidEdgeKind::Crease)
    );
}

#[test]
fn a_group_count_is_checked_against_the_mesh_before_it_is_read() {
    let context = context();
    let triangle = triangle(context);
    let cube = cube(context);
    let ids = [0_u32];
    let mut mesh = entry(triangle, 6, AXIOLID_MESH_ROLE_SHEET);
    mesh.groups = ids.as_ptr();
    mesh.group_mode = AXIOLID_GROUP_MODE_PER_TRIANGLE;
    // The count is far beyond the one readable id, and beyond the mesh.
    for count in [usize::MAX, 12, 2, 0] {
        mesh.group_count = count;
        let (status, result) = Call::new(vec![mesh]).run(context);
        assert_eq!(status, AxiolidStatus::InvalidArgument, "{count}");
        assert_eq!(result, AxiolidResultHandle::INVALID);
        assert_eq!(
            last_error_mesh(context),
            (AxiolidStatus::Ok, 1, 6),
            "{count}"
        );
    }
    // The cube's count on the triangle is refused although the cube accepts it.
    let mut other = entry(cube, 7, AXIOLID_MESH_ROLE_SOLID);
    other.group_mode = AXIOLID_GROUP_MODE_PER_TRIANGLE;
    other.groups = ids.as_ptr();
    other.group_count = 1;
    let (status, _) = Call::new(vec![other]).run(context);
    assert_eq!(status, AxiolidStatus::InvalidArgument);
    assert_eq!(last_error_mesh(context), (AxiolidStatus::Ok, 1, 7));
}

#[test]
fn a_refusal_names_the_mesh_it_concerns() {
    let context = context();
    assert_eq!(last_error_mesh(context).0, AxiolidStatus::NoError);
    let triangle = triangle(context);
    // An open triangle is no solid.
    let call = Call::new(vec![entry(triangle, 9, AXIOLID_MESH_ROLE_SOLID)]);
    let (status, _) = call.run(context);
    assert_eq!(status, AxiolidStatus::OperationFailed);
    assert_eq!(last_error_mesh(context), (AxiolidStatus::Ok, 1, 9));
    let mut info = AxiolidErrorInfo::default();
    assert_eq!(
        unsafe { axiolid_v0_4_context_last_error_info(context, &mut info) },
        AxiolidStatus::Ok
    );
    assert_eq!(info.operation, AxiolidOperation::MeshHiddenLines);
    assert_eq!(info.status, AxiolidStatus::OperationFailed);

    // A refusal of the whole request names no mesh.
    let mut call = Call::new(vec![entry(triangle, 9, AXIOLID_MESH_ROLE_SHEET)]);
    call.limits.max_meshes = 0;
    let (status, _) = call.run(context);
    assert_eq!(status, AxiolidStatus::LimitExceeded);
    assert_eq!(last_error_mesh(context), (AxiolidStatus::Ok, 0, 0));
}

#[test]
fn a_duplicate_id_is_refused_naming_that_id() {
    let context = context();
    let triangle = triangle(context);
    let call = Call::new(vec![
        entry(triangle, 3, AXIOLID_MESH_ROLE_SHEET),
        entry(triangle, 3, AXIOLID_MESH_ROLE_SHEET),
    ]);
    let (status, _) = call.run(context);
    assert_eq!(status, AxiolidStatus::OperationFailed);
    assert_eq!(last_error_mesh(context), (AxiolidStatus::Ok, 1, 3));
}

#[test]
fn invalid_arguments_are_refused_at_the_boundary() {
    let context = context();
    let triangle = triangle(context);
    let base = || Call::new(vec![entry(triangle, 1, AXIOLID_MESH_ROLE_SHEET)]);
    let cases: Vec<(&str, Change)> = vec![
        ("projection", Box::new(|c| c.view.projection = 2)),
        ("view reserved", Box::new(|c| c.view.reserved = [0, 1])),
        ("has_near", Box::new(|c| c.view.has_near = 2)),
        ("has_far", Box::new(|c| c.view.has_far = 2)),
        (
            "origin NaN",
            Box::new(|c| c.view.frame.origin[1] = f64::NAN),
        ),
        (
            "axis infinite",
            Box::new(|c| c.view.frame.z[2] = f64::INFINITY),
        ),
        ("unused near NaN", Box::new(|c| c.view.near = f64::NAN)),
        (
            "unused far infinite",
            Box::new(|c| c.view.far = f64::INFINITY),
        ),
        ("edge reserved", Box::new(|c| c.edges.reserved = [0, 0, 1])),
        ("crease NaN", Box::new(|c| c.edges.crease_cosine = f64::NAN)),
        ("has_crease", Box::new(|c| c.edges.has_crease = 2)),
        ("boundary flag", Box::new(|c| c.edges.boundary = 2)),
        ("silhouette flag", Box::new(|c| c.edges.silhouette = 3)),
        ("hidden flag", Box::new(|c| c.edges.hidden = 255)),
        (
            "group_boundary flag",
            Box::new(|c| c.edges.group_boundary = 2),
        ),
        (
            "region NaN",
            Box::new(|c| {
                c.region = Some(AxiolidImageBox {
                    min: [f64::NAN, 0.0],
                    max: [1.0, 1.0],
                })
            }),
        ),
        ("role", Box::new(|c| c.meshes[0].role = 0)),
        ("group mode", Box::new(|c| c.meshes[0].group_mode = 4)),
        (
            "single with groups",
            Box::new(|c| c.meshes[0].group_count = 1),
        ),
        ("tolerance", Box::new(|c| c.tolerance.linear = f64::NAN)),
    ];
    for (name, change) in cases {
        let mut call = base();
        change(&mut call);
        let (status, result) = call.run(context);
        assert_eq!(status, AxiolidStatus::InvalidArgument, "{name}");
        assert_eq!(result, AxiolidResultHandle::INVALID, "{name}");
        let mut info = AxiolidErrorInfo::default();
        assert_eq!(
            unsafe { axiolid_v0_4_context_last_error_info(context, &mut info) },
            AxiolidStatus::Ok,
            "{name}"
        );
        assert_eq!(info.status, AxiolidStatus::InvalidArgument, "{name}");
        assert_eq!(info.operation, AxiolidOperation::MeshHiddenLines, "{name}");
    }
    let (mut meshes, mut results) = (0_usize, 1_usize);
    unsafe { axiolid_v0_4_context_live_object_counts(context, &mut meshes, &mut results) };
    assert_eq!((meshes, results), (1, 0));
}

#[test]
fn handles_pointers_and_counts_are_checked() {
    let config = AxiolidContextConfig {
        max_meshes: 2,
        ..AxiolidContextConfig::default()
    };
    let context = context_with(config);
    let triangle = triangle(context);
    let call = Call::new(vec![entry(triangle, 1, AXIOLID_MESH_ROLE_SHEET)]);

    let unknown = Call::new(vec![entry(
        AxiolidMeshHandle(u64::MAX),
        5,
        AXIOLID_MESH_ROLE_SHEET,
    )]);
    assert_eq!(unknown.run(context).0, AxiolidStatus::InvalidHandle);
    assert_eq!(last_error_mesh(context), (AxiolidStatus::Ok, 1, 5));
    assert_eq!(
        call.run(AxiolidContextHandle(u64::MAX)).0,
        AxiolidStatus::InvalidHandle
    );

    // Exactly max_meshes entries are allowed.
    let two = Call::new(vec![
        entry(triangle, 1, AXIOLID_MESH_ROLE_SHEET),
        entry(triangle, 2, AXIOLID_MESH_ROLE_SHEET),
    ]);
    assert_eq!(two.run(context).0, AxiolidStatus::Ok);

    // More meshes than the context allows are refused before any is read.
    let many = Call::new(vec![entry(triangle, 1, AXIOLID_MESH_ROLE_SHEET); 3]);
    assert_eq!(many.run(context).0, AxiolidStatus::LimitExceeded);

    unsafe {
        let mut result = AxiolidResultHandle::INVALID;
        let (view, edges, limits) = (plan(), all_edges(), limits());
        let run = |meshes, count, view: *const _, edges: *const _, limits: *const _, out| {
            axiolid_v0_4_hidden_lines(
                context,
                meshes,
                count,
                view,
                edges,
                null(),
                limits,
                TOLERANCE,
                out,
            )
        };
        let entries = [entry(triangle, 1, AXIOLID_MESH_ROLE_SHEET)];
        let ok = |status| assert_eq!(status, AxiolidStatus::NullPointer);
        ok(run(null(), 1, &view, &edges, &limits, &mut result));
        ok(run(
            entries.as_ptr(),
            1,
            null(),
            &edges,
            &limits,
            &mut result,
        ));
        ok(run(
            entries.as_ptr(),
            1,
            &view,
            null(),
            &limits,
            &mut result,
        ));
        ok(run(entries.as_ptr(), 1, &view, &edges, null(), &mut result));
        ok(run(
            entries.as_ptr(),
            1,
            &view,
            &edges,
            &limits,
            std::ptr::null_mut(),
        ));
        // An empty mesh list reaches the contract, which refuses it.
        assert_eq!(
            run(null(), 0, &view, &edges, &limits, &mut result),
            AxiolidStatus::OperationFailed
        );

        // A per-triangle array may not be null.
        let mut bad = entry(triangle, 1, AXIOLID_MESH_ROLE_SHEET);
        bad.group_mode = AXIOLID_GROUP_MODE_PER_TRIANGLE;
        bad.group_count = 1;
        assert_eq!(
            run(&bad, 1, &view, &edges, &limits, &mut result),
            AxiolidStatus::NullPointer
        );

        // Result readers check the pointer, the handle and the kind.
        let (status, drawn) = call.run(context);
        assert_eq!(status, AxiolidStatus::Ok);
        let mut count = 0_usize;
        assert_eq!(
            axiolid_v0_4_hidden_lines_count(context, drawn, std::ptr::null_mut()),
            AxiolidStatus::NullPointer
        );
        assert_eq!(
            axiolid_v0_4_hidden_lines_count(context, AxiolidResultHandle(u64::MAX), &mut count),
            AxiolidStatus::InvalidHandle
        );
        assert_eq!(
            axiolid_v0_4_hidden_lines_copy(context, drawn, std::ptr::null_mut(), count + 1),
            AxiolidStatus::NullPointer
        );
        assert_eq!(
            axiolid_v0_4_hidden_lines_evidence(context, drawn, std::ptr::null_mut()),
            AxiolidStatus::NullPointer
        );
        assert_eq!(
            axiolid_v0_4_context_last_error_mesh(context, std::ptr::null_mut(), &mut 0),
            AxiolidStatus::NullPointer
        );

        let mut exact = AxiolidResultHandle::INVALID;
        assert_eq!(
            axiolid_v0_4_exact_extrude_rectangle(context, 1.0, 1.0, 1.0, TOLERANCE, &mut exact),
            AxiolidStatus::Ok
        );
        assert_eq!(
            axiolid_v0_4_hidden_lines_count(context, exact, &mut count),
            AxiolidStatus::WrongResultKind
        );
    }
}

#[test]
fn the_context_lists_the_hidden_line_capability() {
    unsafe {
        let context = context();
        let mut count = 0_usize;
        assert_eq!(
            axiolid_v0_4_capability_count(context, &mut count),
            AxiolidStatus::Ok
        );
        assert_eq!(count, 6);
        let found = (0..count)
            .map(|index| {
                let mut descriptor = AxiolidCapabilityDescriptor::default();
                assert_eq!(
                    axiolid_v0_4_capability_get(context, index, &mut descriptor),
                    AxiolidStatus::Ok
                );
                descriptor
            })
            .find(|descriptor| descriptor.operation as i32 == 12)
            .expect("hidden lines is listed");
        assert_eq!(found.representation, AxiolidRepresentation::Profile2d);
        let provider = &found.provider[..usize::from(found.provider_len)];
        assert_eq!(provider, b"scalar-hidden-lines");
    }
}
