use super::*;

/// Stable Boolean operand. An integer alias keeps unknown C values defined.
pub type AxiolidBooleanOperator = i32;
pub const AXIOLID_BOOLEAN_UNION: AxiolidBooleanOperator = 1;
pub const AXIOLID_BOOLEAN_INTERSECTION: AxiolidBooleanOperator = 2;
pub const AXIOLID_BOOLEAN_DIFFERENCE: AxiolidBooleanOperator = 3;
pub const AXIOLID_BOOLEAN_SYMMETRIC_DIFFERENCE: AxiolidBooleanOperator = 4;

fn boolean_operator(value: AxiolidBooleanOperator) -> Option<axiolid::core::BooleanOperator> {
    match value {
        AXIOLID_BOOLEAN_UNION => Some(axiolid::core::BooleanOperator::Union),
        AXIOLID_BOOLEAN_INTERSECTION => Some(axiolid::core::BooleanOperator::Intersection),
        AXIOLID_BOOLEAN_DIFFERENCE => Some(axiolid::core::BooleanOperator::Difference),
        AXIOLID_BOOLEAN_SYMMETRIC_DIFFERENCE => {
            Some(axiolid::core::BooleanOperator::SymmetricDifference)
        }
        _ => None,
    }
}

pub(crate) fn insert_result(
    context: &mut Context,
    result: StoredResult,
) -> Result<u64, AxiolidStatus> {
    if context.results.live_count() >= context.config.max_results as usize {
        return Err(AxiolidStatus::LimitExceeded);
    }
    Ok(context.results.insert(result))
}

/// Execute a tolerance-bounded mesh Boolean and return an owned result handle.
/// # Safety
/// Every non-null pointer must be aligned and valid for the documented read or write extent.
#[no_mangle]
pub unsafe extern "C" fn axiolid_v0_4_boolean(
    context: AxiolidContextHandle,
    subject: AxiolidMeshHandle,
    tool: AxiolidMeshHandle,
    operator: AxiolidBooleanOperator,
    tolerance_value: AxiolidTolerance,
    out_result: *mut AxiolidResultHandle,
) -> AxiolidStatus {
    boundary(|| {
        if out_result.is_null() {
            return AxiolidStatus::NullPointer;
        }
        let Ok(tolerance) = tolerance(tolerance_value) else {
            return AxiolidStatus::InvalidArgument;
        };
        let Some(context) = context_entry(context.0) else {
            return AxiolidStatus::InvalidHandle;
        };
        let mut context = lock_unpoisoned(&context);
        let Some(subject) = context.meshes.get(subject.0).cloned() else {
            return AxiolidStatus::InvalidHandle;
        };
        let Some(tool) = context.meshes.get(tool.0).cloned() else {
            return AxiolidStatus::InvalidHandle;
        };
        let Some(operator) = boolean_operator(operator) else {
            return record_error(
                context,
                AxiolidStatus::InvalidArgument,
                AxiolidOperation::MeshBoolean,
                tolerance_value,
                "axiolid-capi",
                "unknown Boolean operator",
            );
        };
        let options = axiolid::contracts::ExecutionOptions::new(tolerance);
        let outcome = match context
            .application
            .boolean(&subject, &tool, operator, &options)
        {
            Ok(outcome) => outcome,
            Err(error) => return record_application_error(context, error),
        };
        let Ok(handle) = insert_result(&mut context, StoredResult::Mesh(outcome.mesh)) else {
            return AxiolidStatus::LimitExceeded;
        };
        // SAFETY: null was rejected; caller provides one writable result handle.
        unsafe { out_result.write(AxiolidResultHandle(handle)) };
        AxiolidStatus::Ok
    })
}

/// Subtract a borrowed array of mesh handles from one subject mesh.
/// # Safety
/// Every non-null pointer must be aligned and valid for the documented read or write extent.
#[no_mangle]
pub unsafe extern "C" fn axiolid_v0_4_subtract_many(
    context: AxiolidContextHandle,
    subject: AxiolidMeshHandle,
    tools: *const AxiolidMeshHandle,
    tool_count: usize,
    tolerance_value: AxiolidTolerance,
    out_result: *mut AxiolidResultHandle,
) -> AxiolidStatus {
    boundary(|| {
        if out_result.is_null() || (tool_count != 0 && tools.is_null()) {
            return AxiolidStatus::NullPointer;
        }
        let Ok(tolerance) = tolerance(tolerance_value) else {
            return AxiolidStatus::InvalidArgument;
        };
        let tool_handles: &[AxiolidMeshHandle] = if tool_count == 0 {
            &[]
        } else {
            // SAFETY: non-zero pointer/length was validated; caller promises a readable array.
            unsafe { std::slice::from_raw_parts(tools, tool_count) }
        };
        let Some(context) = context_entry(context.0) else {
            return AxiolidStatus::InvalidHandle;
        };
        let mut context = lock_unpoisoned(&context);
        let Some(subject) = context.meshes.get(subject.0).cloned() else {
            return AxiolidStatus::InvalidHandle;
        };
        let Some(tools) = tool_handles
            .iter()
            .map(|handle| context.meshes.get(handle.0).cloned())
            .collect::<Option<Vec<_>>>()
        else {
            return AxiolidStatus::InvalidHandle;
        };
        let options = axiolid::contracts::ExecutionOptions::new(tolerance);
        let outcome = match context
            .application
            .subtract_many(&subject, &tools, &options)
        {
            Ok(outcome) => outcome,
            Err(error) => return record_application_error(context, error),
        };
        let Ok(handle) = insert_result(&mut context, StoredResult::Mesh(outcome.mesh)) else {
            return AxiolidStatus::LimitExceeded;
        };
        // SAFETY: null was rejected; caller provides one writable result handle.
        unsafe { out_result.write(AxiolidResultHandle(handle)) };
        AxiolidStatus::Ok
    })
}

/// Construct an exact rectangular prism; no mesh fallback is permitted.
/// # Safety
/// Every non-null pointer must be aligned and valid for the documented read or write extent.
#[no_mangle]
pub unsafe extern "C" fn axiolid_v0_4_exact_extrude_rectangle(
    context: AxiolidContextHandle,
    width: f64,
    height: f64,
    depth: f64,
    tolerance_value: AxiolidTolerance,
    out_result: *mut AxiolidResultHandle,
) -> AxiolidStatus {
    boundary(|| {
        if out_result.is_null() {
            return AxiolidStatus::NullPointer;
        }
        let Ok(tolerance) = tolerance(tolerance_value) else {
            return AxiolidStatus::InvalidArgument;
        };
        if !width.is_finite()
            || !height.is_finite()
            || !depth.is_finite()
            || width <= 0.0
            || height <= 0.0
            || depth <= 0.0
        {
            return AxiolidStatus::InvalidArgument;
        }
        let Some(context) = context_entry(context.0) else {
            return AxiolidStatus::InvalidHandle;
        };
        let mut context = lock_unpoisoned(&context);
        let profile = axiolid::profile::Profile::Rectangle(axiolid::profile::RectangleProfile {
            x: width,
            y: height,
            thickness: None,
            outer_radius: None,
            inner_radius: None,
        });
        let exact = match context.application.extrude_profile_exact(
            &profile,
            axiolid::core::Vec3::Z,
            depth,
            tolerance,
        ) {
            Ok(exact) => exact,
            Err(error) => return record_application_error(context, error),
        };
        let Ok(handle) = insert_result(&mut context, StoredResult::Exact(Box::new(exact))) else {
            return AxiolidStatus::LimitExceeded;
        };
        // SAFETY: null was rejected; caller provides one writable result handle.
        unsafe { out_result.write(AxiolidResultHandle(handle)) };
        AxiolidStatus::Ok
    })
}

/// Exact Boolean is not part of v0.4; this function fails closed by contract.
/// # Safety
/// Every non-null pointer must be aligned and valid for the documented read or write extent.
#[no_mangle]
pub unsafe extern "C" fn axiolid_v0_4_exact_boolean(
    context: AxiolidContextHandle,
    subject: AxiolidMeshHandle,
    tool: AxiolidMeshHandle,
    tolerance_value: AxiolidTolerance,
    out_result: *mut AxiolidResultHandle,
) -> AxiolidStatus {
    boundary(|| {
        if out_result.is_null() {
            return AxiolidStatus::NullPointer;
        }
        if tolerance(tolerance_value).is_err() {
            return AxiolidStatus::InvalidArgument;
        }
        let Some(context) = context_entry(context.0) else {
            return AxiolidStatus::InvalidHandle;
        };
        let context = lock_unpoisoned(&context);
        if context.meshes.get(subject.0).is_none() || context.meshes.get(tool.0).is_none() {
            return AxiolidStatus::InvalidHandle;
        }
        record_error(
            context,
            AxiolidStatus::UnsupportedExact,
            AxiolidOperation::MeshBoolean,
            tolerance_value,
            "none",
            "v0.4 does not provide exact Boolean; no mesh fallback was attempted",
        )
    })
}

/// Query whether a result owns an exact B-rep or triangle mesh.
/// # Safety
/// Every non-null pointer must be aligned and valid for the documented read or write extent.
#[no_mangle]
pub unsafe extern "C" fn axiolid_v0_4_result_kind(
    context: AxiolidContextHandle,
    result: AxiolidResultHandle,
    out_kind: *mut AxiolidGeometryKind,
) -> AxiolidStatus {
    boundary(|| {
        if out_kind.is_null() {
            return AxiolidStatus::NullPointer;
        }
        let Some(context) = context_entry(context.0) else {
            return AxiolidStatus::InvalidHandle;
        };
        let context = lock_unpoisoned(&context);
        let Some(result) = context.results.get(result.0) else {
            return AxiolidStatus::InvalidHandle;
        };
        let kind = match result {
            StoredResult::Mesh(_) => AxiolidGeometryKind::TriangleMesh,
            StoredResult::Exact(exact) => {
                let _ = exact.topology();
                AxiolidGeometryKind::ExactBrep
            }
            StoredResult::HiddenLines(_) => AxiolidGeometryKind::HiddenLines,
            StoredResult::VisibleFaces(_) => AxiolidGeometryKind::VisibleFaces,
        };
        // SAFETY: null was rejected; caller provides one writable enum value.
        unsafe { out_kind.write(kind) };
        AxiolidStatus::Ok
    })
}

/// Consume a mesh result and transfer its mesh to a new context-owned mesh handle.
/// # Safety
/// Every non-null pointer must be aligned and valid for the documented read or write extent.
#[no_mangle]
pub unsafe extern "C" fn axiolid_v0_4_result_take_mesh(
    context: AxiolidContextHandle,
    result: AxiolidResultHandle,
    out_mesh: *mut AxiolidMeshHandle,
) -> AxiolidStatus {
    boundary(|| {
        if out_mesh.is_null() {
            return AxiolidStatus::NullPointer;
        }
        let Some(context) = context_entry(context.0) else {
            return AxiolidStatus::InvalidHandle;
        };
        let mut context = lock_unpoisoned(&context);
        if context.meshes.live_count() >= context.config.max_meshes as usize {
            return AxiolidStatus::LimitExceeded;
        }
        if !matches!(context.results.get(result.0), Some(StoredResult::Mesh(_))) {
            return if context.results.get(result.0).is_some() {
                AxiolidStatus::WrongResultKind
            } else {
                AxiolidStatus::InvalidHandle
            };
        }
        let Some(StoredResult::Mesh(mesh)) = context.results.remove(result.0) else {
            return AxiolidStatus::InvalidHandle;
        };
        let mesh = AxiolidMeshHandle(context.meshes.insert(mesh));
        // SAFETY: null was rejected; caller provides one writable mesh handle.
        unsafe { out_mesh.write(mesh) };
        AxiolidStatus::Ok
    })
}

/// Destroy an operation result and its owned geometry.
#[no_mangle]
pub extern "C" fn axiolid_v0_4_result_destroy(
    context: AxiolidContextHandle,
    result: AxiolidResultHandle,
) -> AxiolidStatus {
    boundary(|| {
        let Some(context) = context_entry(context.0) else {
            return AxiolidStatus::InvalidHandle;
        };
        let mut context = lock_unpoisoned(&context);
        match context.results.remove(result.0) {
            Some(_) => AxiolidStatus::Ok,
            None => AxiolidStatus::InvalidHandle,
        }
    })
}
