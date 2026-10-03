#!/usr/bin/env python3
"""Run copied downstream Rust probes against one immutable Git artifact."""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import NoReturn

import tomllib

ROOT = Path(__file__).resolve().parents[1]
CONSUMERS = ROOT / "tests/consumers"
NATIVE_CONSUMER = ROOT / "tests/native/cmake-consumer"
PROFILES = (
    "linear-intersection-minimal",
    "mesh-rule-checker",
    "2d-curves",
    "parametric-curves",
    "cad-exact",
    "rust-facade-application",
)
REVISION_RE = re.compile(r"[0-9a-f]{40}")
UNVERSIONED_ABI_RE = re.compile(r"\baxiolid_(?!v\d+_\d+_)[A-Za-z0-9_]+")


def fail(message: str) -> NoReturn:
    raise RuntimeError(message)


def run(
    args: list[str],
    *,
    cwd: Path = ROOT,
    env: dict[str, str] | None = None,
    check: bool = True,
) -> subprocess.CompletedProcess[str]:
    # stdout and stderr are captured on separate streams (not merged) so that
    # incidental stderr noise from the invoked process -- e.g. a first-use
    # rustup toolchain-sync notice on a cold CI runner -- can never land in
    # front of machine-readable stdout (JSON from `cargo metadata`) and break
    # a downstream json.loads call. Both streams are still shown on failure.
    process = subprocess.run(
        args,
        cwd=cwd,
        env=env,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if check and process.returncode:
        fail(
            f"{' '.join(args)} failed in {cwd}:\n"
            f"--- stdout ---\n{process.stdout}\n"
            f"--- stderr ---\n{process.stderr}"
        )
    return process


def toml_string(value: str) -> str:
    return json.dumps(value, ensure_ascii=True)


def render_manifest(
    source: str,
    *,
    versions: dict[str, str],
    allowed: set[str],
    source_url: str,
    revision: str,
) -> str:
    """Replace workspace paths with exact-version, exact-revision Git dependencies."""
    if not REVISION_RE.fullmatch(revision):
        raise ValueError("source revision must be 40 lowercase hexadecimal characters")
    if not source_url.startswith(("https://", "file:///")):
        raise ValueError("source URL must be absolute HTTPS or file URI")
    parsed = tomllib.loads(source)
    package = parsed.get("package")
    dependencies = parsed.get("dependencies", {})
    if not isinstance(package, dict) or not isinstance(dependencies, dict):
        raise TypeError("probe manifest must define package and dependencies tables")

    output = ["[package]"]
    for key in ("name", "version", "edition"):
        value = package.get(key)
        if not isinstance(value, str):
            raise TypeError(f"probe package {key} must be a string")
        output.append(f"{key} = {toml_string(value)}")
    output.extend(("publish = false", "", "[workspace]", "", "[dependencies]"))

    for dependency in sorted(dependencies):
        if dependency not in allowed:
            raise ValueError(
                f"probe references unpublished internal package {dependency}"
            )
        version = versions.get(dependency)
        if version is None:
            raise ValueError(f"workspace package version unavailable for {dependency}")
        original = dependencies[dependency]
        if isinstance(original, str):
            original = {}
        if not isinstance(original, dict):
            raise TypeError(f"dependency {dependency} must use a dependency table")
        fields = [
            f"version = {toml_string('=' + version)}",
            f"git = {toml_string(source_url)}",
            f"rev = {toml_string(revision)}",
        ]
        if original.get("default-features") is False:
            fields.append("default-features = false")
        features = original.get("features")
        if features is not None:
            if not isinstance(features, list) or not all(
                isinstance(feature, str) for feature in features
            ):
                raise ValueError(f"dependency {dependency} features must be strings")
            fields.append(
                "features = ["
                + ", ".join(toml_string(value) for value in features)
                + "]"
            )
        output.append(f"{dependency} = {{ {', '.join(fields)} }}")
    output.append("")
    rendered = "\n".join(output)
    if "path" in rendered or "../" in rendered:
        raise ValueError("rendered black-box manifest contains a relative path")
    return rendered


def workspace_packages() -> tuple[dict[str, str], set[str]]:
    metadata = json.loads(
        run(
            ["cargo", "+1.88.0", "metadata", "--no-deps", "--format-version", "1"]
        ).stdout
    )
    versions: dict[str, str] = {}
    allowed: set[str] = set()
    for package in metadata["packages"]:
        manifest = Path(package["manifest_path"])
        try:
            manifest.relative_to(ROOT / "crates")
        except ValueError:
            continue
        name = package["name"]
        versions[name] = package["version"]
        publish = package.get("publish")
        if publish is None or publish:
            allowed.add(name)
    return versions, allowed


def verify_native_source_policy(probe: Path) -> None:
    for path in sorted(probe.rglob("*")):
        if not path.is_file() or path.suffix not in {".c", ".cc", ".cpp", ".h", ".hpp"}:
            continue
        match = UNVERSIONED_ABI_RE.search(path.read_text(encoding="utf-8"))
        if match:
            raise ValueError(f"unversioned ABI symbol {match.group(0)} in {path.name}")


def copy_probe(
    profile: str,
    destination: Path,
    *,
    versions: dict[str, str],
    allowed: set[str],
    source_url: str,
    revision: str,
) -> Path:
    source = CONSUMERS / profile
    probe = destination / profile
    shutil.copytree(
        source, probe, ignore=shutil.ignore_patterns("Cargo.lock", ".gitignore", "target")
    )
    manifest = probe / "Cargo.toml"
    manifest.write_text(
        render_manifest(
            manifest.read_text(encoding="utf-8"),
            versions=versions,
            allowed=allowed,
            source_url=source_url,
            revision=revision,
        ),
        encoding="utf-8",
    )
    return probe


def verify_resolved_boundary(probe: Path, revision: str, env: dict[str, str]) -> None:
    metadata = json.loads(
        run(
            ["cargo", "+1.88.0", "metadata", "--locked", "--format-version", "1"],
            cwd=probe,
            env=env,
        ).stdout
    )
    resolved = [
        package
        for package in metadata["packages"]
        if package["name"].startswith("axiolid") and package.get("source") is not None
    ]
    if not resolved:
        fail(f"{probe.name}: no external Axiolid packages resolved")
    for package in resolved:
        source = package["source"]
        if not source.startswith("git+") or not source.endswith("#" + revision):
            fail(f"{probe.name}: mutable or non-Git Axiolid source {source}")
        if Path(package["manifest_path"]).is_relative_to(ROOT):
            fail(f"{probe.name}: resolved a package from the source workspace")


def run_rust_probes(work: Path, revision: str, env: dict[str, str]) -> None:
    versions, allowed = workspace_packages()
    source_repo = work / "axiolid-release.git"
    run(
        ["git", "clone", "--bare", "--no-local", str(ROOT), str(source_repo)],
        cwd=work,
    )
    run(["git", "cat-file", "-e", revision + "^{commit}"], cwd=source_repo)
    source_url = source_repo.as_uri()
    probes = work / "rust-probes"
    probes.mkdir()
    for profile in PROFILES:
        probe = copy_probe(
            profile,
            probes,
            versions=versions,
            allowed=allowed,
            source_url=source_url,
            revision=revision,
        )
        run(["cargo", "+1.88.0", "generate-lockfile"], cwd=probe, env=env)
        verify_resolved_boundary(probe, revision, env)
        result = run(
            ["cargo", "+1.88.0", "run", "--release", "--locked", "--quiet"],
            cwd=probe,
            env=env,
        )
        if not result.stdout.strip():
            fail(f"{profile}: probe emitted no semantic evidence")
        print(result.stdout.strip())

    facade = probes / "rust-facade-application" / "Cargo.toml"
    mutated = facade.read_text(encoding="utf-8").replace(
        'features = ["application"]', "features = []"
    )
    if mutated == facade.read_text(encoding="utf-8"):
        fail("feature mutation did not alter facade manifest")
    facade.write_text(mutated, encoding="utf-8")
    run(["cargo", "+1.88.0", "generate-lockfile"], cwd=facade.parent, env=env)
    failure = run(
        ["cargo", "+1.88.0", "check", "--locked", "--quiet"],
        cwd=facade.parent,
        env=env,
        check=False,
    )
    if failure.returncode == 0:
        fail("removing the application feature did not break the facade probe")
    print("downstream mutation: required facade feature removal rejected")

    source_manifest = (CONSUMERS / PROFILES[0] / "Cargo.toml").read_text(
        encoding="utf-8"
    )
    first_dependency = next(iter(tomllib.loads(source_manifest)["dependencies"]))
    try:
        render_manifest(
            source_manifest,
            versions=versions,
            allowed=allowed - {first_dependency},
            source_url=source_url,
            revision=revision,
        )
    except ValueError as error:
        if "unpublished" not in str(error):
            raise
    else:
        fail("removing an allowed public package did not break probe rendering")
    print("downstream mutation: unavailable package rejected")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--revision", help="exact source commit; defaults to HEAD")
    args = parser.parse_args()
    try:
        revision = args.revision or run(["git", "rev-parse", "HEAD"]).stdout.strip()
        if not REVISION_RE.fullmatch(revision):
            fail("revision must be a full 40-character lowercase commit")
        verify_native_source_policy(NATIVE_CONSUMER)
        env = os.environ.copy()
        with tempfile.TemporaryDirectory(prefix="axiolid-downstream-") as temporary:
            work = Path(temporary)
            env["CARGO_TARGET_DIR"] = str(
                Path(env.get("CARGO_TARGET_DIR", work / "cargo-target")).resolve()
            )
            run_rust_probes(work, revision, env)
    except (
        OSError,
        RuntimeError,
        TypeError,
        ValueError,
        subprocess.SubprocessError,
    ) as error:
        print(f"downstream consumer probes failed: {error}", file=sys.stderr)
        return 1
    print(f"downstream consumer probes: {revision}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
