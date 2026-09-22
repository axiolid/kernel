#!/usr/bin/env python3
"""Publish Axiolid crates to crates.io in stable dependency order.

`cargo publish --workspace` is unstable in Cargo 1.88. This script derives an
acyclic order from every internal normal, build, and dev dependency, then
publishes one package at a time. Existing immutable versions are skipped only
when their crates.io checksum matches the exact unpatched local package, so a
workflow interrupted by registry propagation or rate limits can be rerun.

The default mode is a side-effect-free plan. Actual publication requires both
`--execute` and Cargo's `CARGO_REGISTRY_TOKEN` environment variable.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import tarfile
import tempfile
import time
import tomllib
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CARGO = os.environ.get("CARGO", "cargo")
USER_AGENT = "axiolid-release-workflow/1.0"


def metadata() -> dict:
    result = subprocess.run(
        [CARGO, "metadata", "--format-version", "1", "--no-deps", "--locked"],
        cwd=ROOT,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    )
    return json.loads(result.stdout)


def publish_plan(data: dict) -> list[dict]:
    members = set(data["workspace_members"])
    packages = {
        package["name"]: package
        for package in data["packages"]
        if package["id"] in members and package.get("publish") != []
    }
    by_path = {
        str(Path(package["manifest_path"]).parent.resolve()): package["name"]
        for package in packages.values()
    }
    dependencies: dict[str, set[str]] = {}
    for name, package in packages.items():
        dependencies[name] = {
            by_path[str(Path(dependency["path"]).resolve())]
            for dependency in package["dependencies"]
            if dependency.get("path")
            and str(Path(dependency["path"]).resolve()) in by_path
        }

    order: list[str] = []
    remaining = {name: set(required) for name, required in dependencies.items()}
    while remaining:
        ready = sorted(name for name, required in remaining.items() if not required)
        if not ready:
            cycle = ", ".join(
                f"{name}->[{', '.join(sorted(required))}]"
                for name, required in sorted(remaining.items())
            )
            raise SystemExit(f"internal publication dependency cycle: {cycle}")
        for name in ready:
            order.append(name)
            del remaining[name]
        published = set(ready)
        for required in remaining.values():
            required.difference_update(published)

    return [packages[name] for name in order]


def registry_checksum(name: str, version: str) -> str | None | bool:
    """Return checksum, False when absent, or None on a transient lookup failure."""
    encoded_name = urllib.parse.quote(name, safe="")
    encoded_version = urllib.parse.quote(version, safe="")
    request = urllib.request.Request(
        f"https://crates.io/api/v1/crates/{encoded_name}/{encoded_version}",
        headers={"User-Agent": USER_AGENT},
    )
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            payload = json.load(response)
            checksum = payload.get("version", {}).get("checksum")
            if not isinstance(checksum, str) or len(checksum) != 64:
                raise SystemExit(f"crates.io returned no valid checksum for {name} {version}")
            return checksum
    except urllib.error.HTTPError as error:
        if error.code == 404:
            return False
        print(f"warning: crates.io lookup for {name} returned HTTP {error.code}", file=sys.stderr)
        return None
    except (OSError, TimeoutError) as error:
        print(f"warning: crates.io lookup for {name} failed: {error}", file=sys.stderr)
        return None


def archive_checksum(archive: Path) -> str:
    digest = hashlib.sha256()
    with archive.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def assert_registry_clean_archive(
    archive: Path, package_name: str, package_version: str, internal_names: set[str]
) -> None:
    with tarfile.open(archive, "r:gz") as package:
        lockfiles = [member for member in package.getmembers() if member.name.endswith("/Cargo.lock")]
        if len(lockfiles) != 1:
            raise SystemExit(
                f"exact upload archive {archive.name} must contain one Cargo.lock, "
                f"found {len(lockfiles)}"
            )
        stream = package.extractfile(lockfiles[0])
        if stream is None:
            raise SystemExit(f"cannot read packaged Cargo.lock from {archive.name}")
        lock = tomllib.loads(stream.read().decode("utf-8"))

    if "patch" in lock:
        raise SystemExit(f"exact upload archive {archive.name} contains patch metadata")
    local_rows = 0
    for row in lock.get("package", []):
        name = row.get("name")
        version = row.get("version")
        source = row.get("source")
        if name == package_name and version == package_version and source is None:
            local_rows += 1
            continue
        if name not in internal_names:
            continue
        checksum = row.get("checksum")
        if source != "registry+https://github.com/rust-lang/crates.io-index":
            raise SystemExit(
                f"exact upload archive {archive.name} resolves internal {name} "
                f"through non-registry source {source!r}"
            )
        if not isinstance(checksum, str) or len(checksum) != 64:
            raise SystemExit(
                f"exact upload archive {archive.name} has no registry checksum for {name}"
            )
    if local_rows != 1:
        raise SystemExit(
            f"exact upload archive {archive.name} has {local_rows} local root lock entries"
        )


def version_already_published(name: str, version: str) -> bool:
    """Whether `name@version` already exists on the registry.

    Distinct from `require_matching_checksum`, which also compares bytes: that
    comparison needs a locally built archive, and building it is exactly the
    cost this check avoids. A published version cannot be replaced, so its
    presence alone justifies skipping.

    `registry_checksum` returns False or None when the crate or version is
    absent, and those cases mean "do the full path" -- never assume published
    on an unclear answer, or a crate would be silently dropped from a release.
    """
    checksum = registry_checksum(name, version)
    return not (checksum is False or checksum is None)


def require_matching_checksum(name: str, version: str, local_checksum: str) -> bool:
    remote_checksum = registry_checksum(name, version)
    if remote_checksum is False or remote_checksum is None:
        return False
    if remote_checksum != local_checksum:
        raise SystemExit(
            f"immutable version collision for {name} {version}: "
            f"local={local_checksum} crates.io={remote_checksum}"
        )
    return True


def wait_for_matching_checksum(
    name: str, version: str, local_checksum: str, args: argparse.Namespace
) -> None:
    for attempt in range(1, args.max_attempts + 1):
        if require_matching_checksum(name, version, local_checksum):
            return
        if attempt == args.max_attempts:
            break
        print(
            f"WAIT_CHECKSUM {name} attempt={attempt + 1} wait={args.dependency_wait}s",
            flush=True,
        )
        time.sleep(args.dependency_wait)
    raise SystemExit(
        f"crates.io did not expose {name} {version} with checksum {local_checksum} "
        f"after {args.max_attempts} attempts"
    )


def retry_delay(output: str, args: argparse.Namespace) -> tuple[int, str] | None:
    lowered = output.lower()
    if "too many requests" in lowered or "429" in lowered:
        return args.rate_limit_wait, "crates.io rate limit"
    if "no matching package named" in lowered or "failed to select a version" in lowered:
        return args.dependency_wait, "registry dependency propagation"
    return None


def reproduce_publish_dry_run(
    package: dict,
    archive: Path,
    environment: dict[str, str],
) -> Path:
    name = package["name"]
    version = package["version"]
    expected_bytes = archive.read_bytes()
    expected_checksum = archive_checksum(archive)
    archive.unlink()
    result = subprocess.run(
        [CARGO, "publish", "-p", name, "--locked", "--dry-run", "--quiet"],
        cwd=ROOT,
        env=environment,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    if result.stdout:
        print(result.stdout, end="", flush=True)
    if result.returncode != 0:
        raise SystemExit(f"cargo publish --dry-run failed for {name} {version}")
    if not archive.is_file():
        raise SystemExit(f"cargo publish --dry-run produced no archive for {name} {version}")
    reproduced_checksum = archive_checksum(archive)
    if reproduced_checksum != expected_checksum or archive.read_bytes() != expected_bytes:
        raise SystemExit(
            f"cargo publish --dry-run archive mismatch for {name} {version}: "
            f"package={expected_checksum} dry-run={reproduced_checksum}"
        )
    print(f"EXACT_ARCHIVE {name} {version} checksum={expected_checksum}", flush=True)
    return archive


def prepare_archive(
    package: dict,
    target: Path,
    environment: dict[str, str],
    args: argparse.Namespace,
) -> Path:
    name = package["name"]
    version = package["version"]
    archive = target / "package" / f"{name}-{version}.crate"
    for attempt in range(1, args.max_attempts + 1):
        archive.unlink(missing_ok=True)
        result = subprocess.run(
            [CARGO, "package", "-p", name, "--locked", "--quiet"],
            cwd=ROOT,
            env=environment,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        if result.stdout:
            print(result.stdout, end="", flush=True)
        if result.returncode == 0:
            if not archive.is_file():
                raise SystemExit(f"cargo package produced no archive for {name} {version}")
            return archive

        retry = retry_delay(result.stdout, args)
        if retry is None:
            raise SystemExit(f"cargo package failed permanently for {name} {version}")
        if attempt == args.max_attempts:
            raise SystemExit(
                f"cargo package exhausted {args.max_attempts} attempts for {name} {version}"
            )
        wait, reason = retry
        print(
            f"RETRY_PACKAGE {name} attempt={attempt + 1} wait={wait}s reason={reason}",
            flush=True,
        )
        time.sleep(wait)
    raise AssertionError("unreachable package retry loop")


def publish(
    package: dict,
    target: Path,
    environment: dict[str, str],
    args: argparse.Namespace,
    internal_names: set[str],
) -> None:
    name = package["name"]
    version = package["version"]
    # Cheap check first. Verifying a crate packages it twice -- once to build
    # the archive, once to prove the publish reproduces the same bytes -- and
    # each packaging compiles the crate. Doing that for every crate on every
    # release is what made a version bump take hours, even for crates whose
    # version was already on the registry and whose upload was then skipped.
    #
    # Registry versions are immutable, so if this exact version exists there is
    # nothing to upload and nothing a byte comparison could change. Skip before
    # the expensive work rather than after it.
    if version_already_published(name, version):
        print(f"SKIP {name} {version}: version already on the registry", flush=True)
        return
    archive = prepare_archive(package, target, environment, args)
    assert_registry_clean_archive(archive, name, version, internal_names)
    archive = reproduce_publish_dry_run(package, archive, environment)
    assert_registry_clean_archive(archive, name, version, internal_names)
    local_checksum = archive_checksum(archive)
    if require_matching_checksum(name, version, local_checksum):
        print(f"SKIP {name} {version}: published checksum matches {local_checksum}", flush=True)
        return

    for attempt in range(1, args.max_attempts + 1):
        result = subprocess.run(
            [CARGO, "publish", "-p", name, "--locked"],
            cwd=ROOT,
            env=environment,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        if result.stdout:
            print(result.stdout, end="", flush=True)
        if result.returncode == 0:
            if not archive.is_file() or archive_checksum(archive) != local_checksum:
                raise SystemExit(
                    f"cargo publish regenerated different local bytes for {name} {version}"
                )
            wait_for_matching_checksum(name, version, local_checksum, args)
            print(f"PUBLISHED {name} {version} checksum={local_checksum}", flush=True)
            return

        output = result.stdout.lower()
        if "already uploaded" in output or "already exists" in output:
            wait_for_matching_checksum(name, version, local_checksum, args)
            print(f"SKIP {name} {version}: published checksum matches", flush=True)
            return

        retry = retry_delay(result.stdout, args)
        if retry is None:
            raise SystemExit(f"cargo publish failed permanently for {name} {version}")
        if attempt == args.max_attempts:
            raise SystemExit(
                f"cargo publish exhausted {args.max_attempts} attempts for {name} {version}"
            )
        wait, reason = retry
        print(
            f"RETRY_PUBLISH {name} attempt={attempt + 1} wait={wait}s reason={reason}",
            flush=True,
        )
        time.sleep(wait)
    raise AssertionError("unreachable publish retry loop")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--execute", action="store_true", help="perform irreversible uploads")
    parser.add_argument("--max-attempts", type=int, default=30)
    parser.add_argument("--dependency-wait", type=int, default=30)
    parser.add_argument("--rate-limit-wait", type=int, default=610)
    args = parser.parse_args()
    if args.max_attempts < 1 or args.dependency_wait < 0 or args.rate_limit_wait < 0:
        parser.error("attempts must be positive and wait values must be non-negative")

    plan = publish_plan(metadata())
    print(f"PUBLISH_PLAN={len(plan)}")
    for index, package in enumerate(plan, 1):
        print(f"{index:02d} {package['name']} {package['version']}")

    if not args.execute:
        print("PLAN_ONLY: pass --execute with CARGO_REGISTRY_TOKEN to publish")
        return
    if not os.environ.get("CARGO_REGISTRY_TOKEN"):
        raise SystemExit("CARGO_REGISTRY_TOKEN is required with --execute")

    with tempfile.TemporaryDirectory(prefix="axiolid-publish-") as temporary:
        target = Path(temporary) / "target"
        environment = os.environ.copy()
        environment["CARGO_TARGET_DIR"] = str(target)
        internal_names = {package["name"] for package in plan}
        for package in plan:
            publish(package, target, environment, args, internal_names)
    print(f"PUBLISH_WORKSPACE=PASS crates={len(plan)}")


if __name__ == "__main__":
    main()
