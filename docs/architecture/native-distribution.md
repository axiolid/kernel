# Native distribution and CMake

Axiolid's native contract is ABI **0.4**. It is independent of Rust's unstable ABI and is exposed only through [`axiolid.h`](https://github.com/axiolid/kernel/blob/main/crates/facade/axiolid-capi/include/axiolid.h). Both consumption paths provide the same CMake target:

```cmake
target_link_libraries(your_target PRIVATE Axiolid::axiolid)
```

Set `AXIOLID_LINKAGE` to `SHARED` (default) or `STATIC` before creating/finding the target.

## Immutable source build

Vendor `native/cmake/AxiolidFetch.cmake`, then pin a full Git object ID—not a branch or mutable tag:

```cmake
include(AxiolidFetch.cmake)
axiolid_fetch(
  GIT_REPOSITORY https://github.com/axiolid/kernel.git
  GIT_COMMIT <40-character-reviewed-commit>
  LINKAGE SHARED
)
```

The helper rejects anything other than a full hexadecimal commit. It builds `axiolid-capi` with Cargo and verifies that the checked-in header matches the Rust API. Rust **1.88.0** and CMake **3.24+** are required. A local `add_subdirectory(<checkout>/native)` path exists for repository development only; it is not the distribution pinning mechanism. The targets set only target-local include/link properties; they never mutate `CMAKE_C_FLAGS`, `CMAKE_CXX_FLAGS`, or inject `target-cpu=native`.

When `native/` is the top-level project, it also supports a conventional install:

```sh
cmake -S native -B build/native -DCMAKE_BUILD_TYPE=Release
cmake --build build/native --config Release
cmake --install build/native --config Release --prefix "$PWD/dist/axiolid"
```

`AXIOLID_ENABLE_INSTALL` defaults off when Axiolid is a subproject, so fetching it does not add files to the consumer's install. Set it explicitly only when that behavior is wanted.

## Release archive

Extract only after verification, then point CMake at its root:

```sh
sha256sum -c axiolid-native-v0.1.1-x86_64-unknown-linux-gnu.tar.gz.sha256
python3 scripts/verify-native-package.py axiolid-native-v0.1.1-x86_64-unknown-linux-gnu.tar.gz --extract-to /opt
cmake -S . -B build -DCMAKE_PREFIX_PATH=/opt/axiolid-native-v0.1.1-x86_64-unknown-linux-gnu
cmake --build build --config Release
```

`find_package(Axiolid 0.1 CONFIG REQUIRED)` then creates:

- `Axiolid::axiolid` — stable selected target;
- `Axiolid::axiolid_shared` — shared implementation;
- `Axiolid::axiolid_static` — static implementation.

On Windows, copy `bin/axiolid_capi.dll` beside the executable when using shared linkage. Unix imported targets retain the extracted library path; deployment may instead install/copy the shared library according to the application's loader policy.

`SHARED` is recommended because it contains Rust's implementation dependencies. `STATIC` is available when the final application accepts system-link responsibility; generated targets add the required Linux/macOS libraries. Windows uses the dynamic UCRT/VCRuntime model (`/MD` family): static Axiolid linkage does not imply `/MT`. v0.4 does not promise a system-wide SONAME—consume the archive through its imported target and use the `axiolid_v0_4_*` prefix as the ABI boundary.

## Supported release targets

| Rust target | Build environment | Native toolchain / runtime | Artifacts | Verification |
|---|---|---|---|---|
| `x86_64-unknown-linux-gnu` | Ubuntu 22.04 | GCC/Clang; glibc 2.35 baseline | `.so`, `.a` | C/C++ Debug/Release execution; shared/static |
| `aarch64-unknown-linux-gnu` | Ubuntu 22.04 x86 cross host | `aarch64-linux-gnu-gcc`; glibc 2.35 sysroot | `.so`, `.a` | cross-build plus ELF/archive machine verification; not executed |
| `x86_64-apple-darwin` | macOS 15 Intel, deployment target 13.0 | Apple Clang; system `libSystem`, Security and CoreFoundation for static linkage | `.dylib`, `.a` | C/C++ Debug/Release execution; shared/static |
| `x86_64-pc-windows-msvc` | Windows Server 2022 | Visual Studio 2022 MSVC; system UCRT/VCRuntime (not a `/MT` bundle) | `.dll`, import `.lib`, static `.lib` | C/C++ Debug/Release execution; shared/static |

For source cross-compilation, set `AXIOLID_CARGO_TARGET` **and** use a matching CMake toolchain file; Cargo's linker must be configured for the same triple. Installed targets encode their OS and processor and fail configuration when a package is selected for the wrong target. Cross compilation never introduces `target-cpu=native`.

Other triples, MinGW, 32-bit processes, and ARM macOS are not claimed by v0.4. The public functions use the platform C calling convention and fixed-width integer fields where width is part of the ABI; `size_t` fields follow the supported 64-bit target ABI. C++ consumers receive `extern "C"` declarations and no C++ runtime dependency.

## Archive contract and integrity

Each `axiolid-native-v0.1.1-<target>` archive contains:

- `include/axiolid.h`;
- shared, static, and platform import libraries;
- `lib/cmake/Axiolid/{AxiolidConfig.cmake,AxiolidConfigVersion.cmake,AxiolidFetch.cmake,AxiolidTargets.cmake}`;
- `LICENSE`;
- `manifest.json` with commit, target, profile, Rust compiler, and per-file SHA-256/size metadata.

Every archive has a `.sha256` sidecar; the release also carries `SHA256SUMS`. `scripts/verify-native-package.py` rejects path traversal, links, missing/extra payload files, dirty release provenance, hash or size drift, wrong binary formats, and wrong object architectures—including members of static archives. Packaging normalizes timestamps, ownership, ordering, and modes; reproducibility is unit-tested and was verified with byte-identical host archives.

## CI and release behavior

`.github/workflows/native.yml` copies out-of-tree C and C++ consumers into a clean temporary tree, builds them through both source and exact extracted-archive paths, and compares their ABI markers and semantic success/typed-refusal behavior. It runs those probes in Debug and Release, with shared and static linkage, on Linux, macOS, and Windows. Mutation runs prove that removing a required symbol, header, package config, or explicit static runtime-link closure breaks the consumer. An additional immutable-source test clones an exact commit through `AxiolidFetch.cmake`. The AArch64 Linux job cross-builds and verifies binary identity but does not claim target execution. Release jobs rebuild exact tagged source, verify every archive, require the complete four-target set, generate consolidated checksums, and upload only after all jobs pass.

Debug and Release select different Cargo profiles but preserve ABI 0.4. Published archives are Release artifacts. A package with an incompatible requested major/minor ABI is rejected by `AxiolidConfigVersion.cmake`.
