#include "axiolid.h"
#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(expr) do { AxiolidStatus s_ = (expr); if (s_ != AxiolidStatus_Ok) { fprintf(stderr, "line %d: status %d\n", __LINE__, (int)s_); return 1; } } while (0)

#define ASSERT_SIZE(type, bytes) _Static_assert(sizeof(type) == (bytes), #type " ABI drift")
ASSERT_SIZE(AxiolidVersion, 12);
ASSERT_SIZE(AxiolidContextHandle, 8);
ASSERT_SIZE(AxiolidMeshHandle, 8);
ASSERT_SIZE(AxiolidResultHandle, 8);
ASSERT_SIZE(AxiolidContextConfig, 20);
ASSERT_SIZE(AxiolidTolerance, 16);
ASSERT_SIZE(AxiolidTransform, 128);
ASSERT_SIZE(AxiolidCapabilityDescriptor, 176);
ASSERT_SIZE(AxiolidBounds, 48);
ASSERT_SIZE(AxiolidMeasurements, 64);
ASSERT_SIZE(AxiolidStatus, 4);
ASSERT_SIZE(AxiolidOperation, 4);
ASSERT_SIZE(AxiolidRepresentation, 4);
ASSERT_SIZE(AxiolidExactness, 4);
ASSERT_SIZE(AxiolidGeometryKind, 4);
ASSERT_SIZE(AxiolidViewFrame, 96);
ASSERT_SIZE(AxiolidView, 120);
ASSERT_SIZE(AxiolidImageBox, 32);
ASSERT_SIZE(AxiolidEdgeSelection, 16);
ASSERT_SIZE(AxiolidEdgeKind, 4);
ASSERT_SIZE(AxiolidVisibility, 4);
ASSERT_SIZE(AxiolidHiddenLineSegment, 56);
#if SIZE_MAX == UINT64_MAX
ASSERT_SIZE(AxiolidViewMesh, 40);
ASSERT_SIZE(AxiolidHiddenLineLimits, 40);
ASSERT_SIZE(AxiolidHiddenLineEvidence, 80);
ASSERT_SIZE(AxiolidMeshAudit, 24);
ASSERT_SIZE(AxiolidErrorInfo, 104);
#endif

static const double CUBE_POSITIONS[] = {
  0,0,0, 2,0,0, 2,2,0, 0,2,0,
  0,0,2, 2,0,2, 2,2,2, 0,2,2
};
static const uint32_t CUBE_INDICES[] = {
  0,2,1, 0,3,2, 4,5,6, 4,6,7, 0,1,5, 0,5,4,
  1,2,6, 1,6,5, 2,3,7, 2,7,6, 3,0,4, 3,4,7
};
static const double TOOL_POSITIONS[] = {
  1,1,0.5, 3,1,0.5, 3,3,0.5, 1,3,0.5,
  1,1,1.5, 3,1,1.5, 3,3,1.5, 1,3,1.5
};

static AxiolidContextConfig config(void) {
  AxiolidContextConfig value = {AXIOLID_PROVIDER_PORTABLE, 64, 64, 1024, 1024};
  return value;
}

static void *context_worker(void *ignored) {
  (void)ignored;
  AxiolidContextHandle context = AxiolidContextHandle_INVALID;
  AxiolidContextConfig limits = config();
  if (axiolid_v0_4_context_create(&limits, &context) != AxiolidStatus_Ok) return (void *)1;
  AxiolidMeshHandle mesh = AxiolidMeshHandle_INVALID;
  if (axiolid_v0_4_mesh_import(context, CUBE_POSITIONS, 8, CUBE_INDICES, 12, &mesh) != AxiolidStatus_Ok) return (void *)1;
  AxiolidTolerance tolerance = {1e-9, 1e-12};
  AxiolidMeasurements measurements = {0};
  if (axiolid_v0_4_mesh_measure(context, mesh, tolerance, &measurements) != AxiolidStatus_Ok) return (void *)1;
  if (fabs(measurements.signed_volume - 8.0) > 1e-9) return (void *)1;
  if (axiolid_v0_4_mesh_destroy(context, mesh) != AxiolidStatus_Ok) return (void *)1;
  if (axiolid_v0_4_context_destroy(context) != AxiolidStatus_Ok) return (void *)1;
  return NULL;
}

int main(void) {
  AxiolidVersion version = {0};
  CHECK(axiolid_v0_4_version(&version));
  if (version.abi_major != 0 || version.abi_minor != 4) return 2;

  AxiolidContextConfig limits = config();
  AxiolidContextHandle context = AxiolidContextHandle_INVALID;
  CHECK(axiolid_v0_4_context_create(&limits, &context));

  size_t capability_count = 0;
  CHECK(axiolid_v0_4_capability_count(context, &capability_count));
  if (capability_count != 5) return 3;

  AxiolidMeshHandle subject = AxiolidMeshHandle_INVALID;
  CHECK(axiolid_v0_4_mesh_import(context, CUBE_POSITIONS, 8, CUBE_INDICES, 12, &subject));
  AxiolidTolerance tolerance = {1e-9, 1e-12};
  AxiolidMeshAudit audit = {0};
  CHECK(axiolid_v0_4_mesh_audit(context, subject, tolerance, &audit));
  if (!audit.closed_two_manifold) return 4;

  AxiolidTransform identity = {{1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1}};
  CHECK(axiolid_v0_4_mesh_transform(context, subject, &identity));
  AxiolidBounds bounds = {0};
  AxiolidMeasurements measurements = {0};
  CHECK(axiolid_v0_4_mesh_bounds(context, subject, &bounds));
  CHECK(axiolid_v0_4_mesh_measure(context, subject, tolerance, &measurements));
  if (fabs(measurements.signed_volume - 8.0) > 1e-9) return 5;

  size_t vertices = 0, triangles = 0;
  CHECK(axiolid_v0_4_mesh_counts(context, subject, &vertices, &triangles));
  double *positions = calloc(vertices * 3, sizeof(double));
  uint32_t *indices = calloc(triangles * 3, sizeof(uint32_t));
  if (!positions || !indices) return 6;
  CHECK(axiolid_v0_4_mesh_copy(context, subject, positions, vertices, indices, triangles));
  free(positions); free(indices);

  AxiolidMeshHandle tool = AxiolidMeshHandle_INVALID;
  CHECK(axiolid_v0_4_mesh_import(context, TOOL_POSITIONS, 8, CUBE_INDICES, 12, &tool));
  AxiolidResultHandle boolean_result = AxiolidResultHandle_INVALID;
  CHECK(axiolid_v0_4_boolean(context, subject, tool, AXIOLID_BOOLEAN_DIFFERENCE, tolerance, &boolean_result));
  AxiolidGeometryKind boolean_kind = AxiolidGeometryKind_Unknown;
  CHECK(axiolid_v0_4_result_kind(context, boolean_result, &boolean_kind));
  if (boolean_kind != AxiolidGeometryKind_TriangleMesh) return 7;
  AxiolidMeshHandle difference = AxiolidMeshHandle_INVALID;
  CHECK(axiolid_v0_4_result_take_mesh(context, boolean_result, &difference));
  AxiolidResultHandle batch_result = AxiolidResultHandle_INVALID;
  CHECK(axiolid_v0_4_subtract_many(context, subject, &tool, 1, tolerance, &batch_result));
  CHECK(axiolid_v0_4_result_destroy(context, batch_result));
  CHECK(axiolid_v0_4_mesh_destroy(context, difference));

  AxiolidResultHandle exact = AxiolidResultHandle_INVALID;
  CHECK(axiolid_v0_4_exact_extrude_rectangle(context, 2, 3, 4, tolerance, &exact));
  AxiolidGeometryKind kind = AxiolidGeometryKind_Unknown;
  CHECK(axiolid_v0_4_result_kind(context, exact, &kind));
  if (kind != AxiolidGeometryKind_ExactBrep) return 7;

  AxiolidResultHandle refused = AxiolidResultHandle_INVALID;
  if (axiolid_v0_4_exact_boolean(context, subject, subject, tolerance, &refused) != AxiolidStatus_UnsupportedExact) return 8;
  AxiolidErrorInfo error = {0};
  CHECK(axiolid_v0_4_context_last_error_info(context, &error));
  if (error.status != AxiolidStatus_UnsupportedExact) return 9;

  double invalid_positions[24];
  memcpy(invalid_positions, CUBE_POSITIONS, sizeof(invalid_positions));
  invalid_positions[0] = NAN;
  AxiolidMeshHandle invalid = AxiolidMeshHandle_INVALID;
  if (axiolid_v0_4_mesh_import(context, invalid_positions, 8, CUBE_INDICES, 12, &invalid) != AxiolidStatus_InvalidArgument) return 10;

  AxiolidViewMesh view_mesh = {subject, 7, NULL, 0, AXIOLID_MESH_ROLE_SOLID, AXIOLID_GROUP_MODE_SINGLE};
  AxiolidView plan = {{{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}}, 0, 0, AXIOLID_PROJECTION_ORTHOGRAPHIC, 0, 0, {0, 0}};
  AxiolidEdgeSelection edges = {0.5, 1, 1, 1, 1, 0, {0, 0, 0}};
  AxiolidHiddenLineLimits line_limits = {8, 1000, 1000, 100000, 1000};
  AxiolidResultHandle drawing = AxiolidResultHandle_INVALID;
  CHECK(axiolid_v0_4_hidden_lines(context, &view_mesh, 1, &plan, &edges, NULL, &line_limits, tolerance, &drawing));
  size_t segment_count = 0;
  CHECK(axiolid_v0_4_hidden_lines_count(context, drawing, &segment_count));
  if (segment_count != 8) return 17;
  AxiolidHiddenLineSegment segments[8];
  if (axiolid_v0_4_hidden_lines_copy(context, drawing, segments, 7) != AxiolidStatus_BufferTooSmall) return 18;
  CHECK(axiolid_v0_4_hidden_lines_copy(context, drawing, segments, 8));
  if (segments[0].mesh_id != 7 || segments[0].kind != AxiolidEdgeKind_Silhouette || segments[0].visibility != AxiolidVisibility_Hidden) return 19;
  AxiolidHiddenLineEvidence evidence = {0};
  CHECK(axiolid_v0_4_hidden_lines_evidence(context, drawing, &evidence));
  if (evidence.visible_segments != 4 || evidence.hidden_segments != 4) return 20;
  view_mesh.id = 9;
  view_mesh.mesh = tool;
  AxiolidResultHandle second_drawing = AxiolidResultHandle_INVALID;
  if (axiolid_v0_4_hidden_lines(context, &view_mesh, 1, &plan, &edges, NULL, &line_limits, tolerance, &second_drawing) != AxiolidStatus_Ok) return 21;
  CHECK(axiolid_v0_4_result_destroy(context, second_drawing));
  plan.projection = 0;
  if (axiolid_v0_4_hidden_lines(context, &view_mesh, 1, &plan, &edges, NULL, &line_limits, tolerance, &second_drawing) != AxiolidStatus_InvalidArgument) return 22;
  uint8_t has_mesh = 7;
  uint64_t named_mesh = 7;
  CHECK(axiolid_v0_4_context_last_error_mesh(context, &has_mesh, &named_mesh));
  if (has_mesh != 0 || named_mesh != 0) return 23;
  CHECK(axiolid_v0_4_result_destroy(context, drawing));

  CHECK(axiolid_v0_4_result_destroy(context, exact));
  if (axiolid_v0_4_result_destroy(context, exact) != AxiolidStatus_InvalidHandle) return 11;
  CHECK(axiolid_v0_4_mesh_destroy(context, tool));
  CHECK(axiolid_v0_4_mesh_destroy(context, subject));
  if (axiolid_v0_4_mesh_destroy(context, subject) != AxiolidStatus_InvalidHandle) return 12;
  size_t live_meshes = 99, live_results = 99;
  CHECK(axiolid_v0_4_context_live_object_counts(context, &live_meshes, &live_results));
  if (live_meshes || live_results) return 13;
  CHECK(axiolid_v0_4_context_destroy(context));
  if (axiolid_v0_4_context_destroy(context) != AxiolidStatus_InvalidHandle) return 14;

  pthread_t threads[4];
  for (size_t i = 0; i < 4; ++i) if (pthread_create(&threads[i], NULL, context_worker, NULL)) return 15;
  for (size_t i = 0; i < 4; ++i) { void *result = NULL; if (pthread_join(threads[i], &result) || result) return 16; }
  puts("axiolid C ABI smoke: ok");
  return 0;
}
