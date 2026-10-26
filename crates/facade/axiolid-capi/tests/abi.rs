use axiolid_capi::{
    axiolid_v0_4_boolean, axiolid_v0_4_capability_count, axiolid_v0_4_capability_get,
    axiolid_v0_4_context_create, axiolid_v0_4_context_destroy,
    axiolid_v0_4_context_last_error_info, axiolid_v0_4_context_last_error_message,
    axiolid_v0_4_context_live_object_counts, axiolid_v0_4_exact_boolean,
    axiolid_v0_4_exact_extrude_rectangle, axiolid_v0_4_mesh_audit, axiolid_v0_4_mesh_bounds,
    axiolid_v0_4_mesh_copy, axiolid_v0_4_mesh_counts, axiolid_v0_4_mesh_destroy,
    axiolid_v0_4_mesh_import, axiolid_v0_4_mesh_measure, axiolid_v0_4_mesh_transform,
    axiolid_v0_4_result_destroy, axiolid_v0_4_result_kind, axiolid_v0_4_result_take_mesh,
    axiolid_v0_4_subtract_many, axiolid_v0_4_version, AxiolidBounds, AxiolidCapabilityDescriptor,
    AxiolidContextConfig, AxiolidContextHandle, AxiolidErrorInfo, AxiolidGeometryKind,
    AxiolidMeasurements, AxiolidMeshAudit, AxiolidMeshHandle, AxiolidResultHandle, AxiolidStatus,
    AxiolidTolerance, AxiolidTransform, AxiolidVersion, AXIOLID_BOOLEAN_DIFFERENCE,
};

#[test]
fn version_and_context_lifecycle_are_stable() {
    unsafe {
        let mut version = AxiolidVersion::default();
        assert_eq!(axiolid_v0_4_version(&mut version), AxiolidStatus::Ok);
        assert_eq!((version.abi_major, version.abi_minor), (0, 4));

        let mut context = AxiolidContextHandle::INVALID;
        let config = AxiolidContextConfig::default();
        assert_eq!(
            axiolid_v0_4_context_create(&config, &mut context),
            AxiolidStatus::Ok
        );
        assert_ne!(context, AxiolidContextHandle::INVALID);
        assert_eq!(axiolid_v0_4_context_destroy(context), AxiolidStatus::Ok);
        assert_eq!(
            axiolid_v0_4_context_destroy(context),
            AxiolidStatus::InvalidHandle
        );
    }
}

fn context() -> AxiolidContextHandle {
    unsafe {
        let mut handle = AxiolidContextHandle::INVALID;
        assert_eq!(
            axiolid_v0_4_context_create(&AxiolidContextConfig::default(), &mut handle),
            AxiolidStatus::Ok
        );
        handle
    }
}

fn cube_arrays() -> (Vec<f64>, Vec<u32>) {
    (
        vec![
            0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 2.0, 2.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0, 2.0, 0.0,
            2.0, 2.0, 2.0, 2.0, 0.0, 2.0, 2.0,
        ],
        vec![
            0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5, 2, 3, 7, 2, 7,
            6, 3, 0, 4, 3, 4, 7,
        ],
    )
}

#[test]
fn capabilities_and_indexed_mesh_round_trip_without_rust_ownership() {
    unsafe {
        let context = context();
        let mut capability_count = 0_usize;
        assert_eq!(
            axiolid_v0_4_capability_count(context, &mut capability_count),
            AxiolidStatus::Ok
        );
        assert_eq!(capability_count, 4);
        let mut descriptor = AxiolidCapabilityDescriptor::default();
        assert_eq!(
            axiolid_v0_4_capability_get(context, 0, &mut descriptor),
            AxiolidStatus::Ok
        );
        assert!(descriptor.id_len > 0);

        let (positions, indices) = cube_arrays();
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
        let (mut vertex_count, mut triangle_count) = (0_usize, 0_usize);
        assert_eq!(
            axiolid_v0_4_mesh_counts(context, mesh, &mut vertex_count, &mut triangle_count),
            AxiolidStatus::Ok
        );
        assert_eq!((vertex_count, triangle_count), (8, 12));

        let mut copied_positions = vec![f64::NAN; vertex_count * 3];
        let mut copied_indices = vec![u32::MAX; triangle_count * 3];
        assert_eq!(
            axiolid_v0_4_mesh_copy(
                context,
                mesh,
                copied_positions.as_mut_ptr(),
                vertex_count,
                copied_indices.as_mut_ptr(),
                triangle_count,
            ),
            AxiolidStatus::Ok
        );
        assert_eq!(copied_positions, positions);
        assert_eq!(copied_indices, indices);

        assert_eq!(axiolid_v0_4_mesh_destroy(context, mesh), AxiolidStatus::Ok);
        assert_eq!(
            axiolid_v0_4_mesh_destroy(context, mesh),
            AxiolidStatus::InvalidHandle
        );
        assert_eq!(axiolid_v0_4_context_destroy(context), AxiolidStatus::Ok);
    }
}

fn import_cube(context: AxiolidContextHandle) -> AxiolidMeshHandle {
    unsafe {
        let (positions, indices) = cube_arrays();
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

#[test]
fn mesh_audit_bounds_measurement_and_transform_are_bounded() {
    unsafe {
        let context = context();
        let mesh = import_cube(context);
        let tolerance = AxiolidTolerance {
            linear: 1.0e-9,
            angular: 1.0e-12,
        };
        let mut audit = AxiolidMeshAudit::default();
        assert_eq!(
            axiolid_v0_4_mesh_audit(context, mesh, tolerance, &mut audit),
            AxiolidStatus::Ok
        );
        assert_eq!(audit.closed_two_manifold, 1);
        assert_eq!(audit.reserved, [0; 7]);

        let mut bounds = AxiolidBounds::default();
        assert_eq!(
            axiolid_v0_4_mesh_bounds(context, mesh, &mut bounds),
            AxiolidStatus::Ok
        );
        assert_eq!(bounds.min, [0.0; 3]);
        assert_eq!(bounds.max, [2.0; 3]);

        let mut measurements = AxiolidMeasurements::default();
        assert_eq!(
            axiolid_v0_4_mesh_measure(context, mesh, tolerance, &mut measurements),
            AxiolidStatus::Ok
        );
        assert!((measurements.surface_area - 24.0).abs() < 1.0e-10);
        assert!((measurements.signed_volume - 8.0).abs() < 1.0e-10);

        let transform = AxiolidTransform {
            columns: [
                1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 5.0, -2.0, 3.0, 1.0,
            ],
        };
        assert_eq!(
            axiolid_v0_4_mesh_transform(context, mesh, &transform),
            AxiolidStatus::Ok
        );
        assert_eq!(
            axiolid_v0_4_mesh_bounds(context, mesh, &mut bounds),
            AxiolidStatus::Ok
        );
        assert_eq!(bounds.min, [5.0, -2.0, 3.0]);
        assert_eq!(bounds.max, [7.0, 0.0, 5.0]);

        assert_eq!(axiolid_v0_4_mesh_destroy(context, mesh), AxiolidStatus::Ok);
        assert_eq!(axiolid_v0_4_context_destroy(context), AxiolidStatus::Ok);
    }
}

fn import_shifted_cube(context: AxiolidContextHandle) -> AxiolidMeshHandle {
    unsafe {
        let (mut positions, indices) = cube_arrays();
        for point in positions.as_chunks_mut::<3>().0 {
            point[0] += 1.0;
            point[1] += 1.0;
            point[2] = point[2] * 0.5 + 0.5;
        }
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

#[test]
fn operation_results_have_explicit_kind_and_transfer_ownership() {
    unsafe {
        let context = context();
        let subject = import_cube(context);
        let tool = import_shifted_cube(context);
        let tolerance = AxiolidTolerance {
            linear: 1.0e-9,
            angular: 1.0e-12,
        };
        let mut result = AxiolidResultHandle::INVALID;
        assert_eq!(
            axiolid_v0_4_boolean(context, subject, tool, 999, tolerance, &mut result),
            AxiolidStatus::InvalidArgument
        );
        assert_eq!(
            axiolid_v0_4_boolean(
                context,
                subject,
                tool,
                AXIOLID_BOOLEAN_DIFFERENCE,
                tolerance,
                &mut result,
            ),
            AxiolidStatus::Ok
        );
        let mut kind = AxiolidGeometryKind::Unknown;
        assert_eq!(
            axiolid_v0_4_result_kind(context, result, &mut kind),
            AxiolidStatus::Ok
        );
        assert_eq!(kind, AxiolidGeometryKind::TriangleMesh);
        let mut output_mesh = AxiolidMeshHandle::INVALID;
        assert_eq!(
            axiolid_v0_4_result_take_mesh(context, result, &mut output_mesh),
            AxiolidStatus::Ok
        );
        assert_eq!(
            axiolid_v0_4_result_destroy(context, result),
            AxiolidStatus::InvalidHandle
        );
        assert_eq!(
            axiolid_v0_4_mesh_destroy(context, output_mesh),
            AxiolidStatus::Ok
        );

        let tools = [tool];
        assert_eq!(
            axiolid_v0_4_subtract_many(
                context,
                subject,
                tools.as_ptr(),
                tools.len(),
                tolerance,
                &mut result,
            ),
            AxiolidStatus::Ok
        );
        assert_eq!(
            axiolid_v0_4_result_destroy(context, result),
            AxiolidStatus::Ok
        );

        assert_eq!(
            axiolid_v0_4_exact_extrude_rectangle(context, 2.0, 3.0, 4.0, tolerance, &mut result),
            AxiolidStatus::Ok
        );
        assert_eq!(
            axiolid_v0_4_result_kind(context, result, &mut kind),
            AxiolidStatus::Ok
        );
        assert_eq!(kind, AxiolidGeometryKind::ExactBrep);
        let mut wrong_kind_mesh = AxiolidMeshHandle::INVALID;
        assert_eq!(
            axiolid_v0_4_result_take_mesh(context, result, &mut wrong_kind_mesh),
            AxiolidStatus::WrongResultKind
        );
        assert_eq!(wrong_kind_mesh, AxiolidMeshHandle::INVALID);
        assert_eq!(
            axiolid_v0_4_result_destroy(context, result),
            AxiolidStatus::Ok
        );

        assert_eq!(
            axiolid_v0_4_exact_boolean(context, subject, tool, tolerance, &mut result),
            AxiolidStatus::UnsupportedExact
        );
        assert_eq!(
            axiolid_v0_4_mesh_destroy(context, subject),
            AxiolidStatus::Ok
        );
        assert_eq!(axiolid_v0_4_mesh_destroy(context, tool), AxiolidStatus::Ok);
        assert_eq!(axiolid_v0_4_context_destroy(context), AxiolidStatus::Ok);
    }
}

#[test]
fn refusal_and_invalid_input_publish_context_owned_structured_errors() {
    unsafe {
        let context = context();
        let subject = import_cube(context);
        let tool = import_shifted_cube(context);
        let tolerance = AxiolidTolerance {
            linear: 1.0e-9,
            angular: 1.0e-12,
        };
        let mut result = AxiolidResultHandle::INVALID;
        assert_eq!(
            axiolid_v0_4_exact_boolean(context, subject, tool, tolerance, &mut result),
            AxiolidStatus::UnsupportedExact
        );
        let mut info = AxiolidErrorInfo::default();
        assert_eq!(
            axiolid_v0_4_context_last_error_info(context, &mut info),
            AxiolidStatus::Ok
        );
        assert_eq!(info.status, AxiolidStatus::UnsupportedExact);
        assert_eq!(info.operation as i32, 6);
        assert_eq!(info.tolerance, tolerance);
        assert_eq!(info.reserved, [0; 7]);

        let mut required = 0_usize;
        assert_eq!(
            axiolid_v0_4_context_last_error_message(
                context,
                std::ptr::null_mut(),
                0,
                &mut required
            ),
            AxiolidStatus::BufferTooSmall
        );
        let mut message = vec![0_u8; required];
        assert_eq!(
            axiolid_v0_4_context_last_error_message(
                context,
                message.as_mut_ptr(),
                message.len(),
                &mut required,
            ),
            AxiolidStatus::Ok
        );
        assert!(std::str::from_utf8(&message[..message.len() - 1])
            .unwrap()
            .contains("no mesh fallback"));

        let positions = [f64::NAN, 0.0, 0.0];
        let indices: [u32; 0] = [];
        let mut mesh = AxiolidMeshHandle::INVALID;
        assert_eq!(
            axiolid_v0_4_mesh_import(
                context,
                positions.as_ptr(),
                1,
                indices.as_ptr(),
                0,
                &mut mesh,
            ),
            AxiolidStatus::InvalidArgument
        );
        assert_eq!(
            axiolid_v0_4_context_last_error_info(context, &mut info),
            AxiolidStatus::Ok
        );
        assert_eq!(info.status, AxiolidStatus::InvalidArgument);

        assert_eq!(
            axiolid_v0_4_mesh_destroy(context, subject),
            AxiolidStatus::Ok
        );
        assert_eq!(axiolid_v0_4_mesh_destroy(context, tool), AxiolidStatus::Ok);
        assert_eq!(axiolid_v0_4_context_destroy(context), AxiolidStatus::Ok);
    }
}

#[test]
fn live_counts_prove_cleanup_and_independent_contexts_are_thread_safe() {
    unsafe {
        let first = context();
        let second = context();
        let mesh = import_cube(first);
        let (mut vertices, mut triangles) = (0_usize, 0_usize);
        assert_eq!(
            axiolid_v0_4_mesh_counts(second, mesh, &mut vertices, &mut triangles),
            AxiolidStatus::InvalidHandle
        );
        assert_eq!(axiolid_v0_4_mesh_destroy(first, mesh), AxiolidStatus::Ok);
        assert_eq!(axiolid_v0_4_context_destroy(first), AxiolidStatus::Ok);
        assert_eq!(axiolid_v0_4_context_destroy(second), AxiolidStatus::Ok);

        let workers: Vec<_> =
            (0..8)
                .map(|_| {
                    std::thread::spawn(|| {
                        let context = context();
                        let mesh = import_cube(context);
                        let (mut meshes, mut results) = (0_usize, 0_usize);
                        assert_eq!(
                            axiolid_v0_4_context_live_object_counts(
                                context,
                                &mut meshes,
                                &mut results,
                            ),
                            AxiolidStatus::Ok
                        );
                        assert_eq!((meshes, results), (1, 0));
                        assert_eq!(axiolid_v0_4_mesh_destroy(context, mesh), AxiolidStatus::Ok);
                        assert_eq!(
                            axiolid_v0_4_context_live_object_counts(
                                context,
                                &mut meshes,
                                &mut results,
                            ),
                            AxiolidStatus::Ok
                        );
                        assert_eq!((meshes, results), (0, 0));
                        assert_eq!(axiolid_v0_4_context_destroy(context), AxiolidStatus::Ok);
                    })
                })
                .collect();
        for worker in workers {
            worker.join().expect("ABI worker must not panic");
        }
    }
}

#[test]
fn context_budgets_and_index_validation_fail_before_ownership_changes() {
    unsafe {
        let mut invalid_context = AxiolidContextHandle::INVALID;
        let invalid_config = AxiolidContextConfig {
            max_vertices_per_mesh: 0,
            ..AxiolidContextConfig::default()
        };
        assert_eq!(
            axiolid_v0_4_context_create(&invalid_config, &mut invalid_context),
            AxiolidStatus::InvalidArgument
        );
        let unknown_provider = AxiolidContextConfig {
            provider_profile: 999,
            ..AxiolidContextConfig::default()
        };
        assert_eq!(
            axiolid_v0_4_context_create(&unknown_provider, &mut invalid_context),
            AxiolidStatus::InvalidArgument
        );

        let mut context_handle = AxiolidContextHandle::INVALID;
        let config = AxiolidContextConfig {
            provider_profile: axiolid_capi::AXIOLID_PROVIDER_PORTABLE,
            max_vertices_per_mesh: 8,
            max_triangles_per_mesh: 12,
            max_meshes: 3,
            max_results: 1,
        };
        assert_eq!(
            axiolid_v0_4_context_create(&config, &mut context_handle),
            AxiolidStatus::Ok
        );
        let mut empty = AxiolidMeshHandle::INVALID;
        assert_eq!(
            axiolid_v0_4_mesh_import(
                context_handle,
                std::ptr::null(),
                0,
                std::ptr::null(),
                0,
                &mut empty,
            ),
            AxiolidStatus::Ok
        );
        assert_eq!(
            axiolid_v0_4_mesh_copy(
                context_handle,
                empty,
                std::ptr::null_mut(),
                0,
                std::ptr::null_mut(),
                0,
            ),
            AxiolidStatus::Ok
        );
        assert_eq!(
            axiolid_v0_4_mesh_destroy(context_handle, empty),
            AxiolidStatus::Ok
        );
        let subject = import_cube(context_handle);
        let tool = import_shifted_cube(context_handle);
        let (positions, indices) = cube_arrays();
        let invalid_indices = [0_u32, 1, 99];
        let mut rejected = AxiolidMeshHandle::INVALID;
        assert_eq!(
            axiolid_v0_4_mesh_import(
                context_handle,
                positions.as_ptr(),
                positions.len() / 3,
                invalid_indices.as_ptr(),
                1,
                &mut rejected,
            ),
            AxiolidStatus::InvalidArgument
        );
        assert_eq!(rejected, AxiolidMeshHandle::INVALID);
        let third = import_cube(context_handle);
        assert_eq!(
            axiolid_v0_4_mesh_import(
                context_handle,
                positions.as_ptr(),
                positions.len() / 3,
                indices.as_ptr(),
                indices.len() / 3,
                &mut rejected,
            ),
            AxiolidStatus::LimitExceeded
        );

        let tolerance = AxiolidTolerance {
            linear: 1.0e-9,
            angular: 1.0e-12,
        };
        let mut first = AxiolidResultHandle::INVALID;
        let mut second = AxiolidResultHandle::INVALID;
        assert_eq!(
            axiolid_v0_4_subtract_many(
                context_handle,
                subject,
                std::ptr::null(),
                0,
                tolerance,
                &mut first,
            ),
            AxiolidStatus::Ok
        );
        assert_ne!(first, AxiolidResultHandle::INVALID);
        assert_eq!(
            axiolid_v0_4_result_destroy(context_handle, first),
            AxiolidStatus::Ok
        );
        first = AxiolidResultHandle::INVALID;
        assert_eq!(
            axiolid_v0_4_boolean(
                context_handle,
                subject,
                tool,
                AXIOLID_BOOLEAN_DIFFERENCE,
                tolerance,
                &mut first,
            ),
            AxiolidStatus::Ok
        );
        assert_eq!(
            axiolid_v0_4_boolean(
                context_handle,
                subject,
                tool,
                AXIOLID_BOOLEAN_DIFFERENCE,
                tolerance,
                &mut second,
            ),
            AxiolidStatus::LimitExceeded
        );
        assert_eq!(second, AxiolidResultHandle::INVALID);
        assert_eq!(
            axiolid_v0_4_result_destroy(context_handle, first),
            AxiolidStatus::Ok
        );
        assert_eq!(
            axiolid_v0_4_mesh_destroy(context_handle, third),
            AxiolidStatus::Ok
        );
        assert_eq!(
            axiolid_v0_4_mesh_destroy(context_handle, subject),
            AxiolidStatus::Ok
        );
        assert_eq!(
            axiolid_v0_4_mesh_destroy(context_handle, tool),
            AxiolidStatus::Ok
        );
        assert_eq!(
            axiolid_v0_4_context_destroy(context_handle),
            AxiolidStatus::Ok
        );
    }
}
