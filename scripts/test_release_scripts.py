#!/usr/bin/env python3
"""Release-script regression tests; network-free and side-effect-free."""

from __future__ import annotations

import hashlib
import importlib.util
import io
from pathlib import Path
import re
import tarfile
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]


def load_script(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / filename)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {filename}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


publish = load_script("publish_workspace", "publish-workspace.py")
verify = load_script("verify_packages", "verify-packages.py")
prepare_release = load_script("prepare_release", "prepare-release.py")


def write_crate_archive(path: Path, lock_text: str | None) -> None:
    with tarfile.open(path, "w:gz") as archive:
        payloads = {"crate-0.1.0/src/lib.rs": b"pub fn probe() {}\n"}
        if lock_text is not None:
            payloads["crate-0.1.0/Cargo.lock"] = lock_text.encode()
        for name, payload in payloads.items():
            member = tarfile.TarInfo(name)
            member.size = len(payload)
            member.mtime = 0
            archive.addfile(member, io.BytesIO(payload))


class PublishWorkspaceTests(unittest.TestCase):
    def test_workflow_actions_are_pinned_to_immutable_commits(self) -> None:
        workflows = sorted((ROOT / ".github" / "workflows").glob("*.y*ml"))
        self.assertTrue(workflows)
        action_refs = []
        for workflow in workflows:
            action_refs.extend(
                (workflow, match.group(1))
                for match in re.finditer(r"uses:\s+[^@\s]+@([^#\s]+)", workflow.read_text())
            )
        self.assertTrue(action_refs)
        for workflow, action_ref in action_refs:
            self.assertRegex(
                action_ref,
                r"^[0-9a-f]{40}$",
                f"{workflow.relative_to(ROOT)} uses mutable action ref {action_ref}",
            )

    def test_plan_orders_every_internal_dependency_before_its_consumer(self) -> None:
        data = publish.metadata()
        plan = publish.publish_plan(data)
        positions = {package["name"]: index for index, package in enumerate(plan)}
        by_path = {
            str(Path(package["manifest_path"]).parent.resolve()): package["name"]
            for package in plan
        }
        self.assertEqual(len(plan), 50)
        for package in plan:
            for dependency in package["dependencies"]:
                dependency_path = dependency.get("path")
                if dependency_path and str(Path(dependency_path).resolve()) in by_path:
                    dependency_name = by_path[str(Path(dependency_path).resolve())]
                    self.assertLess(
                        positions[dependency_name],
                        positions[package["name"]],
                        f"{dependency_name} must publish before {package['name']}",
                    )

    def test_matching_registry_checksum_is_idempotent(self) -> None:
        checksum = "a" * 64
        with mock.patch.object(publish, "registry_checksum", return_value=checksum):
            self.assertTrue(publish.require_matching_checksum("crate", "1.0.0", checksum))

    def test_version_collision_fails_closed(self) -> None:
        with mock.patch.object(publish, "registry_checksum", return_value="b" * 64):
            with self.assertRaisesRegex(SystemExit, "immutable version collision"):
                publish.require_matching_checksum("crate", "1.0.0", "a" * 64)

    def test_absent_registry_version_is_not_treated_as_published(self) -> None:
        with mock.patch.object(publish, "registry_checksum", return_value=False):
            self.assertFalse(publish.require_matching_checksum("crate", "1.0.0", "a" * 64))

    def test_registry_checksum_treats_http_429_as_transient(self) -> None:
        error = publish.urllib.error.HTTPError(
            "https://crates.io", 429, "Too Many Requests", None, None
        )
        with mock.patch.object(publish.urllib.request, "urlopen", side_effect=error):
            with mock.patch("sys.stderr") as stderr:
                self.assertIsNone(publish.registry_checksum("crate", "1.0.0"))
        self.assertTrue(stderr.write.called)

    def test_registry_checksum_treats_transport_error_as_transient(self) -> None:
        with mock.patch.object(
            publish.urllib.request, "urlopen", side_effect=OSError("network unavailable")
        ):
            with mock.patch("sys.stderr") as stderr:
                self.assertIsNone(publish.registry_checksum("crate", "1.0.0"))
        self.assertTrue(stderr.write.called)

    def test_archive_checksum_hashes_exact_bytes(self) -> None:
        payload = b"normalized crate archive\x00"
        with tempfile.TemporaryDirectory() as temporary:
            archive = Path(temporary) / "crate-1.0.0.crate"
            archive.write_bytes(payload)
            self.assertEqual(publish.archive_checksum(archive), hashlib.sha256(payload).hexdigest())

    def test_prepare_archive_uses_unpatched_locked_package(self) -> None:
        package = {"name": "leaf", "version": "0.1.0"}
        args = SimpleNamespace(max_attempts=1, dependency_wait=0, rate_limit_wait=0)
        with tempfile.TemporaryDirectory() as temporary:
            target = Path(temporary) / "target"

            def fake_run(command, **kwargs):
                self.assertEqual(
                    command,
                    [publish.CARGO, "package", "-p", "leaf", "--locked", "--quiet"],
                )
                self.assertNotIn("--config", command)
                self.assertEqual(kwargs["env"]["CARGO_TARGET_DIR"], str(target))
                archive = target / "package" / "leaf-0.1.0.crate"
                archive.parent.mkdir(parents=True)
                archive.write_bytes(b"crate")
                return SimpleNamespace(returncode=0, stdout="")

            environment = {"CARGO_TARGET_DIR": str(target)}
            with mock.patch.object(publish.subprocess, "run", side_effect=fake_run):
                archive = publish.prepare_archive(package, target, environment, args)
            self.assertEqual(archive.read_bytes(), b"crate")

    def test_bootstrap_archive_rejects_packaged_lock_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            archive = Path(temporary) / "crate-0.1.0.crate"
            write_crate_archive(archive, "version = 4\n")
            with self.assertRaisesRegex(SystemExit, "contains Cargo.lock/patch metadata"):
                verify.assert_no_packaged_lockfile(archive)

    def test_exact_archive_rejects_patch_metadata(self) -> None:
        lock = """version = 4
[[package]]
name = "crate"
version = "0.1.0"
[[patch.unused]]
name = "internal"
version = "0.1.0"
"""
        with tempfile.TemporaryDirectory() as temporary:
            archive = Path(temporary) / "crate-0.1.0.crate"
            write_crate_archive(archive, lock)
            with self.assertRaisesRegex(SystemExit, "contains patch metadata"):
                publish.assert_registry_clean_archive(
                    archive, "crate", "0.1.0", {"crate", "internal"}
                )

    def test_exact_archive_rejects_path_backed_internal_dependency(self) -> None:
        lock = """version = 4
[[package]]
name = "crate"
version = "0.1.0"
[[package]]
name = "internal"
version = "0.1.0"
"""
        with tempfile.TemporaryDirectory() as temporary:
            archive = Path(temporary) / "crate-0.1.0.crate"
            write_crate_archive(archive, lock)
            with self.assertRaisesRegex(SystemExit, "non-registry source"):
                publish.assert_registry_clean_archive(
                    archive, "crate", "0.1.0", {"crate", "internal"}
                )

    def test_publish_dry_run_must_reproduce_archive_bytes(self) -> None:
        package = {"name": "crate", "version": "0.1.0"}
        with tempfile.TemporaryDirectory() as temporary:
            target = Path(temporary) / "target" / "package"
            target.mkdir(parents=True)
            archive = target / "crate-0.1.0.crate"
            archive.write_bytes(b"package bytes")

            def fake_run(command, **kwargs):
                self.assertEqual(
                    command,
                    [publish.CARGO, "publish", "-p", "crate", "--locked", "--dry-run", "--quiet"],
                )
                archive.write_bytes(b"different dry-run bytes")
                return SimpleNamespace(returncode=0, stdout="")

            with mock.patch.object(publish.subprocess, "run", side_effect=fake_run):
                with self.assertRaisesRegex(SystemExit, "archive mismatch"):
                    publish.reproduce_publish_dry_run(package, archive, {})

    def test_prepare_archive_fails_closed_on_permanent_error(self) -> None:
        package = {"name": "leaf", "version": "0.1.0"}
        args = SimpleNamespace(max_attempts=1, dependency_wait=0, rate_limit_wait=0)
        result = SimpleNamespace(returncode=101, stdout="")
        with tempfile.TemporaryDirectory() as temporary:
            target = Path(temporary) / "target"
            with mock.patch.object(publish.subprocess, "run", return_value=result):
                with self.assertRaisesRegex(SystemExit, "failed permanently"):
                    publish.prepare_archive(package, target, {}, args)


    def test_prepare_release_requires_a_forward_version_bump(self) -> None:
        with self.assertRaisesRegex(SystemExit, "strictly greater"):
            prepare_release.require_forward_bump("0.1.0", "0.1.0")
        with self.assertRaisesRegex(SystemExit, "strictly greater"):
            prepare_release.require_forward_bump("0.2.0", "0.1.9")
        prepare_release.require_forward_bump("0.1.0", "0.2.0")

    def test_prepare_release_rejects_malformed_semver(self) -> None:
        with self.assertRaisesRegex(SystemExit, "not a valid semantic version"):
            prepare_release.parse_semver("v0.2.0")
        with self.assertRaisesRegex(SystemExit, "not a valid semantic version"):
            prepare_release.parse_semver("0.2")

    def test_prepare_release_rolls_unreleased_into_a_dated_heading(self) -> None:
        prefix = "# Changelog\n\n"
        body = "\n### Added\n- one\n\n"
        suffix = "## [0.1.0] - 2026-01-01\n\nfirst\n"
        rolled = prepare_release.rolled_changelog(prefix, body, suffix, "0.2.0", "2026-09-03")
        self.assertIn("## [Unreleased]\n\n## [0.2.0] - 2026-09-03\n### Added\n- one\n\n", rolled)
        self.assertTrue(rolled.endswith(suffix))

    def test_prepare_release_rejects_an_empty_unreleased_section(self) -> None:
        with self.assertRaisesRegex(SystemExit, "nothing to release"):
            prepare_release.require_nonempty_unreleased("\nno entries here\n")

    def test_prepare_release_bumps_workspace_and_every_internal_dependency(self) -> None:
        original = prepare_release.CARGO_TOML
        try:
            with tempfile.TemporaryDirectory() as temporary:
                cargo_toml = Path(temporary) / "Cargo.toml"
                cargo_toml.write_text(
                    "[workspace.package]\n"
                    'version = "0.1.0"\n'
                    "\n"
                    "[workspace.dependencies]\n"
                    'axiolid-core = { path = "crates/foundation/core", version = "0.1.0" }\n'
                    'glam = "0.29"\n',
                    encoding="utf-8",
                )
                prepare_release.CARGO_TOML = cargo_toml
                bumped = prepare_release.bumped_workspace_toml("0.2.0", "0.1.0")
                self.assertIn('version = "0.2.0"\n', bumped)
                self.assertIn(
                    'axiolid-core = { path = "crates/foundation/core", version = "0.2.0" }',
                    bumped,
                )
                self.assertIn('glam = "0.29"', bumped)
        finally:
            prepare_release.CARGO_TOML = original

    def test_prepare_release_fails_closed_when_no_dependency_matches_current_version(
        self,
    ) -> None:
        original = prepare_release.CARGO_TOML
        try:
            with tempfile.TemporaryDirectory() as temporary:
                cargo_toml = Path(temporary) / "Cargo.toml"
                cargo_toml.write_text(
                    '[workspace.package]\nversion = "0.1.0"\n\n[workspace.dependencies]\nglam = "0.29"\n',
                    encoding="utf-8",
                )
                prepare_release.CARGO_TOML = cargo_toml
                with self.assertRaisesRegex(SystemExit, "no internal axiolid-\\* dependency"):
                    prepare_release.bumped_workspace_toml("0.2.0", "0.1.0")
        finally:
            prepare_release.CARGO_TOML = original


if __name__ == "__main__":
    unittest.main()
