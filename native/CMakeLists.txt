cmake_minimum_required(VERSION 3.24)
project(AxiolidNative VERSION 0.4.0 LANGUAGES C)

if(NOT CMAKE_SIZEOF_VOID_P EQUAL 8)
  message(FATAL_ERROR "Axiolid ABI v0.4 supports only 64-bit consumers")
endif()
if(WIN32 AND NOT MSVC)
  message(FATAL_ERROR "Axiolid ABI v0.4 Windows builds require MSVC; MinGW is not supported")
endif()
if(NOT WIN32 AND NOT APPLE AND NOT CMAKE_SYSTEM_NAME STREQUAL "Linux")
  message(FATAL_ERROR "Axiolid ABI v0.4 supports Linux, macOS, and Windows only")
endif()
set(_install_default OFF)
if(PROJECT_IS_TOP_LEVEL)
  set(_install_default ON)
endif()
option(AXIOLID_ENABLE_INSTALL "Install the Axiolid native package" ${_install_default})

set(AXIOLID_SOURCE_ROOT "${CMAKE_CURRENT_LIST_DIR}/.." CACHE PATH "Axiolid source root")
set(AXIOLID_CARGO_TARGET "" CACHE STRING "Optional Rust target triple")
set(AXIOLID_LINKAGE "SHARED" CACHE STRING "Library behind Axiolid::axiolid (SHARED or STATIC)")
set_property(CACHE AXIOLID_LINKAGE PROPERTY STRINGS SHARED STATIC)
string(TOUPPER "${AXIOLID_LINKAGE}" AXIOLID_LINKAGE)
if(NOT AXIOLID_LINKAGE MATCHES "^(SHARED|STATIC)$")
  message(FATAL_ERROR "AXIOLID_LINKAGE must be SHARED or STATIC")
endif()

find_program(AXIOLID_CARGO_EXECUTABLE cargo REQUIRED)
set(AXIOLID_CARGO_TARGET_DIR "${CMAKE_CURRENT_BINARY_DIR}/cargo-target" CACHE PATH "Cargo target directory")
set(_target_args)
set(_artifact_root "${AXIOLID_CARGO_TARGET_DIR}")
if(AXIOLID_CARGO_TARGET)
  list(APPEND _target_args --target "${AXIOLID_CARGO_TARGET}")
  string(APPEND _artifact_root "/${AXIOLID_CARGO_TARGET}")
endif()

if(WIN32)
  set(_shared_name "axiolid_capi.dll")
  set(_shared_implib_name "axiolid_capi.dll.lib")
  set(_static_name "axiolid_capi.lib")
  set(_axiolid_implib_byproducts
    "${_artifact_root}/debug/${_shared_implib_name}"
    "${_artifact_root}/release/${_shared_implib_name}"
  )
elseif(APPLE)
  set(_shared_name "libaxiolid_capi.dylib")
  set(_static_name "libaxiolid_capi.a")
  set(_axiolid_implib_byproducts)
else()
  set(_shared_name "libaxiolid_capi.so")
  set(_static_name "libaxiolid_capi.a")
  set(_axiolid_implib_byproducts)
endif()

if(PROJECT_IS_TOP_LEVEL AND NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

add_custom_target(axiolid-cargo-build ALL
  COMMAND "${CMAKE_COMMAND}" -E env "CARGO_TARGET_DIR=${AXIOLID_CARGO_TARGET_DIR}"
          "${AXIOLID_CARGO_EXECUTABLE}" xtask ffi check
  COMMAND "${CMAKE_COMMAND}" -E env "CARGO_TARGET_DIR=${AXIOLID_CARGO_TARGET_DIR}"
          "${AXIOLID_CARGO_EXECUTABLE}" build
          $<$<NOT:$<CONFIG:Debug>>:--release>
          --manifest-path "${AXIOLID_SOURCE_ROOT}/Cargo.toml"
          --package axiolid-capi ${_target_args}
  WORKING_DIRECTORY "${AXIOLID_SOURCE_ROOT}"
  BYPRODUCTS
    "${_artifact_root}/debug/${_shared_name}"
    "${_artifact_root}/release/${_shared_name}"
    "${_artifact_root}/debug/${_static_name}"
    "${_artifact_root}/release/${_static_name}"
    ${_axiolid_implib_byproducts}
  COMMAND_EXPAND_LISTS
  VERBATIM
)

add_library(axiolid_shared SHARED IMPORTED GLOBAL)
add_library(Axiolid::axiolid_shared ALIAS axiolid_shared)
add_dependencies(axiolid_shared axiolid-cargo-build)
set_target_properties(axiolid_shared PROPERTIES
  IMPORTED_LOCATION_DEBUG "${_artifact_root}/debug/${_shared_name}"
  IMPORTED_LOCATION_RELEASE "${_artifact_root}/release/${_shared_name}"
  IMPORTED_LOCATION_RELWITHDEBINFO "${_artifact_root}/release/${_shared_name}"
  IMPORTED_LOCATION_MINSIZEREL "${_artifact_root}/release/${_shared_name}"
  INTERFACE_INCLUDE_DIRECTORIES "${AXIOLID_SOURCE_ROOT}/crates/facade/axiolid-capi/include"
)
if(WIN32)
  set_target_properties(axiolid_shared PROPERTIES
    IMPORTED_IMPLIB_DEBUG "${_artifact_root}/debug/${_shared_implib_name}"
    IMPORTED_IMPLIB_RELEASE "${_artifact_root}/release/${_shared_implib_name}"
    IMPORTED_IMPLIB_RELWITHDEBINFO "${_artifact_root}/release/${_shared_implib_name}"
    IMPORTED_IMPLIB_MINSIZEREL "${_artifact_root}/release/${_shared_implib_name}"
  )
elseif(APPLE)
  # The dylib's install name is `@rpath/libaxiolid_capi.dylib` (set in
  # axiolid-capi/build.rs so apps can relocate it). dyld resolves `@rpath`
  # only through the consumer's LC_RPATH (#113). Observed in CI (otool -l):
  # a consumer of the INSTALLED package gets an LC_RPATH to its lib/ from
  # CMake without help, but a consumer of this source-tree target gets none
  # and fails at load. So only this target supplies the dylib's directory.
  set_property(TARGET axiolid_shared APPEND PROPERTY INTERFACE_LINK_OPTIONS
    "LINKER:-rpath,${_artifact_root}/$<IF:$<CONFIG:Debug>,debug,release>")
endif()

add_library(axiolid_static STATIC IMPORTED GLOBAL)
add_library(Axiolid::axiolid_static ALIAS axiolid_static)
add_dependencies(axiolid_static axiolid-cargo-build)
set_target_properties(axiolid_static PROPERTIES
  IMPORTED_LOCATION_DEBUG "${_artifact_root}/debug/${_static_name}"
  IMPORTED_LOCATION_RELEASE "${_artifact_root}/release/${_static_name}"
  IMPORTED_LOCATION_RELWITHDEBINFO "${_artifact_root}/release/${_static_name}"
  IMPORTED_LOCATION_MINSIZEREL "${_artifact_root}/release/${_static_name}"
  INTERFACE_INCLUDE_DIRECTORIES "${AXIOLID_SOURCE_ROOT}/crates/facade/axiolid-capi/include"
)
if(UNIX AND NOT APPLE)
  set_property(TARGET axiolid_static APPEND PROPERTY INTERFACE_LINK_LIBRARIES "dl;pthread;m")
elseif(APPLE)
  set_property(TARGET axiolid_static APPEND PROPERTY INTERFACE_LINK_LIBRARIES "-framework Security;-framework CoreFoundation")
elseif(WIN32)
  # rustc --print=native-static-libs for axiolid-capi on x86_64-pc-windows-msvc
  # (cross-compiled, dev and release profiles, confirmed identical): std pulls
  # these in for env/threading/RNG/error-reporting. msvcrt is handled by
  # CMake's own runtime-library selection, not listed here.
  set_property(TARGET axiolid_static APPEND PROPERTY INTERFACE_LINK_LIBRARIES "kernel32.lib;ntdll.lib;userenv.lib;ws2_32.lib;dbghelp.lib")
endif()

add_library(axiolid_native INTERFACE)
add_library(Axiolid::axiolid ALIAS axiolid_native)
if(AXIOLID_LINKAGE STREQUAL "STATIC")
  target_link_libraries(axiolid_native INTERFACE Axiolid::axiolid_static)
else()
  target_link_libraries(axiolid_native INTERFACE Axiolid::axiolid_shared)
endif()

if(AXIOLID_ENABLE_INSTALL)
  if(WIN32)
    set(SHARED_LOCATION "bin/${_shared_name}")
    set(SHARED_IMPLIB_PROPERTY "set_target_properties(Axiolid::axiolid_shared PROPERTIES IMPORTED_IMPLIB \"\${PACKAGE_PREFIX_DIR}/lib/${_shared_implib_name}\")")
    set(_shared_destination bin)
    set(_implib "${_artifact_root}/$<IF:$<CONFIG:Debug>,debug,release>/${_shared_implib_name}")
  else()
    set(SHARED_LOCATION "lib/${_shared_name}")
    set(SHARED_IMPLIB_PROPERTY "")
    set(_shared_destination lib)
    set(_implib "")
  endif()
  set(STATIC_LOCATION "lib/${_static_name}")
  set(SYSTEM_NAME "${CMAKE_SYSTEM_NAME}")
  if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64|arm64|ARM64)$")
    set(PROCESSOR_REGEX "^(aarch64|arm64|ARM64)$")
  else()
    set(PROCESSOR_REGEX "^(x86_64|AMD64|amd64)$")
  endif()
  if(UNIX AND NOT APPLE)
    set(STATIC_LINK_LIBRARIES "dl;pthread;m")
  elseif(APPLE)
    set(STATIC_LINK_LIBRARIES "-framework CoreFoundation;-framework Security")
  else()
    # rustc --print=native-static-libs for axiolid-capi on
    # x86_64-pc-windows-msvc (cross-compiled, dev and release profiles,
    # confirmed identical): std pulls these in for env/threading/RNG/
    # error-reporting. msvcrt is handled by CMake's own runtime-library
    # selection, not listed here.
    set(STATIC_LINK_LIBRARIES "kernel32.lib;ntdll.lib;userenv.lib;ws2_32.lib;dbghelp.lib")
  endif()
  configure_file(
    "${CMAKE_CURRENT_LIST_DIR}/cmake/AxiolidTargets.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/AxiolidTargets.cmake" @ONLY)
  set(_profile "$<IF:$<CONFIG:Debug>,debug,release>")
  install(FILES
    "${AXIOLID_SOURCE_ROOT}/crates/facade/axiolid-capi/include/axiolid.h"
    DESTINATION include)
  install(FILES "${_artifact_root}/${_profile}/${_shared_name}"
    DESTINATION "${_shared_destination}")
  install(FILES "${_artifact_root}/${_profile}/${_static_name}"
    DESTINATION lib)
  if(_implib)
    install(FILES "${_implib}" DESTINATION lib)
  endif()
  install(FILES
    "${CMAKE_CURRENT_LIST_DIR}/cmake/AxiolidConfig.cmake"
    "${CMAKE_CURRENT_LIST_DIR}/cmake/AxiolidConfigVersion.cmake"
    "${CMAKE_CURRENT_LIST_DIR}/cmake/AxiolidFetch.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/AxiolidTargets.cmake"
    DESTINATION lib/cmake/Axiolid)
  install(FILES "${AXIOLID_SOURCE_ROOT}/LICENSE" DESTINATION .)
endif()
