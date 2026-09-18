cmake_minimum_required(VERSION 3.24)
project(AxiolidCMakeConsumer LANGUAGES C CXX)

if(AXIOLID_SOURCE_TREE)
  add_subdirectory("${AXIOLID_SOURCE_TREE}/native" axiolid-native)
elseif(AXIOLID_PACKAGE_ROOT)
  find_package(Axiolid 0.4 CONFIG REQUIRED
    PATHS "${AXIOLID_PACKAGE_ROOT}"
    NO_DEFAULT_PATH)
else()
  find_package(Axiolid 0.4 CONFIG REQUIRED)
endif()

add_executable(axiolid-cmake-consumer main.c)
add_executable(axiolid-cmake-consumer-cpp main.cpp)
foreach(consumer IN ITEMS axiolid-cmake-consumer axiolid-cmake-consumer-cpp)
  target_link_libraries(${consumer} PRIVATE Axiolid::axiolid)

  # A consumer that calls an undeclared C function must FAIL, not warn.
  # MSVC reports C4013 and GCC/Clang before 14 report a warning, so the
  # symbol-removal mutation could still link on those toolchains and the
  # gate passed while proving nothing (kernel#56).
  if(MSVC)
    target_compile_options(${consumer} PRIVATE $<$<COMPILE_LANGUAGE:C>:/we4013>)
  else()
    target_compile_options(${consumer} PRIVATE
      $<$<COMPILE_LANGUAGE:C>:-Werror=implicit-function-declaration>)
  endif()

  # And the link must resolve eagerly. Mach-O ld and MSVC can otherwise
  # defer or tolerate a missing definition that GNU ld rejects outright,
  # which is the other half of the same blind spot.
  if(APPLE)
    target_link_options(${consumer} PRIVATE LINKER:-undefined,error)
  endif()
  if(WIN32 AND AXIOLID_LINKAGE STREQUAL "SHARED")
    add_custom_command(TARGET ${consumer} POST_BUILD
      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
              "$<TARGET_FILE:Axiolid::axiolid_shared>"
              "$<TARGET_FILE_DIR:${consumer}>"
    )
  endif()
endforeach()

enable_testing()
add_test(NAME axiolid-cmake-consumer COMMAND axiolid-cmake-consumer)
add_test(NAME axiolid-cmake-consumer-cpp COMMAND axiolid-cmake-consumer-cpp)
