#!/usr/bin/env python3
"""Prove source-build and packaged CMake consumers have equivalent behavior."""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CONSUMER = ROOT / "tests/native/cmake-consumer"
C_MARKER = re.compile(r"axiolid native consumer: (\d+\.\d+)")
CPP_MARKER = re.compile(r"axiolid native C\+\+ consumer: (\d+\.\d+)")


def run(args: list[str], *, env: dict[str, str] | None = None) -> str:
    process = subprocess.run(
        args,
        cwd=ROOT,
        env=env,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    if process.returncode:
        print(process.stdout, file=sys.stderr)
        raise RuntimeError(f"command failed ({process.returncode}): {' '.join(args)}")
    return process.stdout


def configure_build_test(
    consumer: Path,
    build: Path,
    build_type: str,
    linkage: str,
    extra: list[str],
    env: dict[str, str],
) -> str:
    run(
        [
            "cmake",
            "-S",
            str(consumer),
            "-B",
            str(build),
            f"-DCMAKE_BUILD_TYPE={build_type}",
            f"-DAXIOLID_LINKAGE={linkage}",
            *extra,
        ],
        env=env,
    )
    run(
        ["cmake", "--build", str(build), "--config", build_type, "--parallel", "2"],
        env=env,
    )
    output = run(
        [
            "ctest",
            "--test-dir",
            str(build),
            "-C",
            build_type,
            "-V",
            "--output-on-failure",
        ],
        env=env,
    )
    c_match = C_MARKER.search(output)
    cpp_match = CPP_MARKER.search(output)
    if not c_match or not cpp_match:
        raise RuntimeError("native C and C++ consumers did not both emit ABI markers")
    if c_match.group(1) != cpp_match.group(1):
        raise RuntimeError("native C and C++ consumers reported different ABI versions")
    return c_match.group(1)


def assert_embedded_build_type_unchanged(
    work: Path, cargo_target: Path, env: dict[str, str]
) -> None:
    source = work / "embedding"
    source.mkdir()
    (source / "CMakeLists.txt").write_text("""cmake_minimum_required(VERSION 3.24)
project(Embedding C)
set(before "${CMAKE_BUILD_TYPE}")
add_subdirectory("${AXIOLID_ROOT}/native" axiolid)
if(NOT "${CMAKE_BUILD_TYPE}" STREQUAL "${before}")
  message(FATAL_ERROR "Axiolid changed the consumer CMAKE_BUILD_TYPE")
endif()
""")
    run(
        [
            "cmake",
            "-S",
            str(source),
            "-B",
            str(work / "embedding-build"),
            f"-DAXIOLID_ROOT={ROOT}",
            f"-DAXIOLID_CARGO_TARGET_DIR={cargo_target}",
        ],
        env=env,
    )


def install_source_package(
    build: Path,
    prefix: Path,
    build_type: str,
    cargo_target: Path,
    env: dict[str, str],
) -> None:
    run(
        [
            "cmake",
            "-S",
            str(ROOT / "native"),
            "-B",
            str(build),
            f"-DCMAKE_BUILD_TYPE={build_type}",
            f"-DAXIOLID_CARGO_TARGET_DIR={cargo_target}",
        ],
        env=env,
    )
    run(
        ["cmake", "--build", str(build), "--config", build_type, "--parallel", "2"],
        env=env,
    )
    run(
        [
            "cmake",
            "--install",
            str(build),
            "--config",
            build_type,
            "--prefix",
            str(prefix),
        ],
        env=env,
    )


def expect_consumer_failure(
    consumer: Path,
    build: Path,
    build_type: str,
    linkage: str,
    package_root: Path,
    env: dict[str, str],
    *,
    configure_may_fail: bool,
) -> None:
    configured = subprocess.run(
        [
            "cmake",
            "-S",
            str(consumer),
            "-B",
            str(build),
            f"-DCMAKE_BUILD_TYPE={build_type}",
            f"-DAXIOLID_LINKAGE={linkage}",
            f"-DAXIOLID_PACKAGE_ROOT={package_root}",
        ],
        env=env,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    if configured.returncode:
        if configure_may_fail:
            return
        raise RuntimeError(
            "symbol mutation failed during configuration instead of compile/link"
        )
    built = subprocess.run(
        ["cmake", "--build", str(build), "--config", build_type, "--parallel", "2"],
        env=env,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    if built.returncode == 0:
        raise RuntimeError(f"mutated consumer unexpectedly built: {consumer.name}")


def assert_strict_undeclared_call(consumer: Path) -> None:
    """The symbol mutation only proves anything if calling an undeclared
    function is a hard error.

    MSVC reports C4013 and older GCC/Clang report a warning, so without an
    explicit promotion the mutated consumer still compiles and links and the
    gate passes while proving nothing (kernel#56). Assert the flags are
    present rather than trusting that every CI toolchain errors by default.
    """
    text = (consumer / "CMakeLists.txt").read_text(encoding="utf-8")
    for needed in ("/we4013", "-Werror=implicit-function-declaration"):
        if needed not in text:
            raise RuntimeError(
                f"consumer must promote implicit declarations to errors: {needed} missing"
            )
    # The Mach-O eager-resolution flag is asserted for the same reason, and it
    # needs asserting MORE than the others: no Linux run exercises it, so
    # deleting it leaves every local and ubuntu CI job green while the macOS
    # gate quietly stops proving anything. That is the exact shape of the
    # original kernel#56 blind spot, so it is checked by text rather than by
    # a platform that happens to be running.
    if "-undefined,error" not in text:
        raise RuntimeError(
            "consumer must ask Mach-O to resolve undefined symbols eagerly: "
            "-undefined,error missing (no non-Apple job can catch this)"
        )


def assert_package_mutations(
    consumer: Path,
    package_root: Path,
    work: Path,
    build_type: str,
    linkage: str,
    env: dict[str, str],
) -> None:
    symbol_consumer = work / "mutation-symbol-consumer"
    shutil.copytree(consumer, symbol_consumer)
    for name in ("main.c", "main.cpp"):
        source = symbol_consumer / name
        mutated = source.read_text(encoding="utf-8").replace(
            "axiolid_v0_4_version", "axiolid_v0_4_removed_symbol"
        )
        if mutated == source.read_text(encoding="utf-8"):
            raise RuntimeError(f"symbol mutation did not alter {name}")
        source.write_text(mutated, encoding="utf-8")
    expect_consumer_failure(
        symbol_consumer,
        work / "mutation-symbol-build",
        build_type,
        linkage,
        package_root,
        env,
        configure_may_fail=False,
    )
    no_header = work / "mutation-no-header-package"
    shutil.copytree(package_root, no_header)
    (no_header / "include/axiolid.h").unlink()
    expect_consumer_failure(
        consumer,
        work / "mutation-header-build",
        build_type,
        linkage,
        no_header,
        env,
        configure_may_fail=True,
    )

    no_config = work / "mutation-no-config-package"
    shutil.copytree(package_root, no_config)
    (no_config / "lib/cmake/Axiolid/AxiolidConfig.cmake").unlink()
    expect_consumer_failure(
        consumer,
        work / "mutation-package-build",
        build_type,
        linkage,
        no_config,
        env,
        configure_may_fail=True,
    )

    mutations = ["symbol", "header", "package"]
    if linkage == "STATIC":
        no_runtime_links = work / "mutation-no-runtime-links-package"
        shutil.copytree(package_root, no_runtime_links)
        targets = no_runtime_links / "lib/cmake/Axiolid/AxiolidTargets.cmake"
        original = targets.read_text(encoding="utf-8")
        mutated = re.sub(
            r'^  INTERFACE_LINK_LIBRARIES "[^\"]+"\n',
            "",
            original,
            count=1,
            flags=re.MULTILINE,
        )
        if mutated == original and os.name != "nt":
            raise RuntimeError(
                "static package declares no runtime-link closure to test"
            )
        if mutated != original:
            targets.write_text(mutated, encoding="utf-8")
            expect_consumer_failure(
                consumer,
                work / "mutation-runtime-links-build",
                build_type,
                linkage,
                no_runtime_links,
                env,
                configure_may_fail=False,
            )
            mutations.append("runtime links")
    print(f"native package mutations rejected: {', '.join(mutations)}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-type", choices=("Debug", "Release"), default="Release")
    parser.add_argument("--linkage", choices=("SHARED", "STATIC"), default="SHARED")
    parser.add_argument(
        "--mutations",
        action="store_true",
        help="prove required symbol, header, and package removal are rejected",
    )
    args = parser.parse_args()
    profile = "debug" if args.build_type == "Debug" else "release"
    try:
        with tempfile.TemporaryDirectory(prefix="axiolid-cmake-") as temporary:
            work = Path(temporary)
            env = os.environ.copy()
            cargo_target = Path(
                env.get("CARGO_TARGET_DIR", work / "cargo-target")
            ).resolve()
            env["CARGO_TARGET_DIR"] = str(cargo_target)
            consumer = work / "consumer"
            shutil.copytree(CONSUMER, consumer)
            assert_embedded_build_type_unchanged(work, cargo_target, env)
            source = configure_build_test(
                consumer,
                work / "source",
                args.build_type,
                args.linkage,
                [
                    f"-DAXIOLID_SOURCE_TREE={ROOT}",
                    f"-DAXIOLID_CARGO_TARGET_DIR={cargo_target}",
                ],
                env,
            )
            install_prefix = work / "installed"
            install_source_package(
                work / "install-build",
                install_prefix,
                args.build_type,
                cargo_target,
                env,
            )
            installed = configure_build_test(
                consumer,
                work / "installed-consumer",
                args.build_type,
                args.linkage,
                [f"-DAXIOLID_PACKAGE_ROOT={install_prefix}"],
                env,
            )
            package_output = work / "packages"
            package_log = run(
                [
                    sys.executable,
                    str(ROOT / "scripts/package-native.py"),
                    "--profile",
                    profile,
                    "--output-dir",
                    str(package_output),
                    "--allow-dirty",
                ],
                env=env,
            )
            archive = Path(package_log.strip().splitlines()[-1])
            extract_to = work / "extracted"
            verify_log = run(
                [
                    sys.executable,
                    str(ROOT / "scripts/verify-native-package.py"),
                    str(archive),
                    "--allow-dirty",
                    "--extract-to",
                    str(extract_to),
                ],
                env=env,
            )
            package_root = Path(verify_log.strip().splitlines()[-1])
            packaged = configure_build_test(
                consumer,
                work / "package",
                args.build_type,
                args.linkage,
                [f"-DAXIOLID_PACKAGE_ROOT={package_root}"],
                env,
            )
            if source != installed or source != packaged or source != "0.4":
                raise RuntimeError(
                    "source/install/archive behavior differs: "
                    f"{source!r}, {installed!r}, {packaged!r}"
                )
            if args.mutations:
                assert_strict_undeclared_call(consumer)
                assert_package_mutations(
                    consumer,
                    package_root,
                    work,
                    args.build_type,
                    args.linkage,
                    env,
                )
    except (OSError, RuntimeError, subprocess.SubprocessError) as error:
        print(f"native CMake test failed: {error}", file=sys.stderr)
        return 1
    print(f"native CMake equivalence: {args.build_type} {args.linkage} ABI {source}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
