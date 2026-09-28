# Native integration

CMake support for `axiolid-capi`, the C ABI. A C or C++ project links the one
stable target `Axiolid::axiolid`, in any of three ways:

- a source build through `add_subdirectory(<kernel>/native)` (this
  `CMakeLists.txt`, which drives Cargo),
- an installed or extracted release archive through
  `find_package(Axiolid 0.4 CONFIG)` (`cmake/AxiolidConfig.cmake`),
- `cmake/AxiolidFetch.cmake`, which fetches a source tree pinned to an
  immutable 40-hex `GIT_COMMIT` and refuses a branch or tag.

`AXIOLID_LINKAGE` selects `SHARED` (default) or `STATIC`. The consumer guide is
[downstream integration](../docs/guide/downstream-integration.md). The ABI is
specified in [C ABI v0.4](../docs/architecture/c-abi-v0.4.md), and archives
in [native distribution](../docs/architecture/native-distribution.md).

## Design notes

- Platform library file names (`.so`, `.dylib`, `.dll` and import library) are
  set in one block of `CMakeLists.txt`. Keep them there.
- Nothing here may set host-specific codegen or the consumer's global
  C/C++ flags. `scripts/check-native-packaging.sh` rejects both, and it
  scans this directory, so name the flags only in that script.
- Platform-specific CI setup belongs in `.github/workflows/native.yml`. The
  CMake files and test fixtures stay platform-neutral.

## Tests

`tests/native/cmake-consumer/` is a black-box C and C++ consumer. It uses only
`find_package` or `add_subdirectory`, the `Axiolid::axiolid` target and the
generated public header, and it must run both a successful operation and a
typed-refusal path. `scripts/test-native-cmake.py` copies it outside the
workspace and builds it from the source tree, from an installed package and
from an extracted archive. With `--mutations`, it also shows that a missing
required symbol, header or package config breaks the build.

```bash
scripts/check-native-packaging.sh    # everything below, as the gate runs it
python3 scripts/test-native-cmake.py --build-type Release --linkage SHARED --mutations
python3 -m unittest tests/native/test_native_packaging.py
```

`tests/native/test_native_packaging.py` covers the archive layouts, fail-closed
archive paths, reproducible writers and binary checks. `tests/native/reject-mutable-ref.cmake` shows that
`axiolid_fetch` refuses a mutable ref.
