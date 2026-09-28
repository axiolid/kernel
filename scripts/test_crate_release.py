#!/usr/bin/env python3
"""Regression tests for prepare-crate-release.py.

Network-free and side-effect-free: no subprocess, no cargo metadata call, only
the pure functions the script exposes, mirroring test_release_scripts.py's
pattern for prepare-release.py. The per-crate changelog page it feeds is
assembled and tested by `cargo xtask docs` (tools/xtask/src/docs/changelog.rs).
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


def load_script(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / filename)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {filename}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


prepare_crate = load_script("prepare_crate_release", "prepare-crate-release.py")


class PrepareCrateReleaseTests(unittest.TestCase):
    def test_requires_a_forward_version_bump(self) -> None:
        with self.assertRaisesRegex(SystemExit, "strictly greater"):
            prepare_crate.require_forward_bump("0.3.0", "0.3.0")
        with self.assertRaisesRegex(SystemExit, "strictly greater"):
            prepare_crate.require_forward_bump("0.3.1", "0.3.0")
        prepare_crate.require_forward_bump("0.3.0", "0.3.1")

    def test_patch_bump_on_0x_is_compatible(self) -> None:
        # 0.x: the minor is the breaking slot, so a patch bump is always fine.
        prepare_crate.require_compatible_bump("0.3.0", "0.3.1")
        prepare_crate.require_compatible_bump("0.3.4", "0.3.99")

    def test_minor_bump_on_0x_is_refused_as_breaking(self) -> None:
        with self.assertRaisesRegex(SystemExit, "breaking bump"):
            prepare_crate.require_compatible_bump("0.3.0", "0.4.0")

    def test_major_bump_on_0x_is_refused_as_breaking(self) -> None:
        with self.assertRaisesRegex(SystemExit, "breaking bump"):
            prepare_crate.require_compatible_bump("0.3.0", "1.0.0")

    def test_patch_bump_past_1_0_is_compatible(self) -> None:
        # 1.x+: the major is the breaking slot, minor/patch are both fine.
        prepare_crate.require_compatible_bump("1.2.3", "1.2.4")
        prepare_crate.require_compatible_bump("1.2.3", "1.9.0")

    def test_major_bump_past_1_0_is_refused_as_breaking(self) -> None:
        with self.assertRaisesRegex(SystemExit, "breaking bump"):
            prepare_crate.require_compatible_bump("1.2.3", "2.0.0")

    def test_bumped_crate_toml_replaces_only_the_version_line(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            manifest = Path(temporary) / "Cargo.toml"
            manifest.write_text(
                "[package]\n"
                'name = "axiolid-probe"\n'
                'version = "0.3.0"\n'
                "\n"
                "[dependencies]\n"
                'axiolid-core = { path = "../core", version = "0.3" }\n',
                encoding="utf-8",
            )
            bumped = prepare_crate.bumped_crate_toml(manifest, "0.3.0", "0.3.1")
            self.assertIn('version = "0.3.1"\n', bumped)
            # A dependency requirement that happens to read "0.3" must survive
            # untouched: only the package's OWN version line may change.
            self.assertIn('version = "0.3" }', bumped)

    def test_bumped_crate_toml_fails_closed_on_version_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            manifest = Path(temporary) / "Cargo.toml"
            manifest.write_text('version = "0.9.0"\n', encoding="utf-8")
            with self.assertRaisesRegex(SystemExit, "did not contain the expected"):
                prepare_crate.bumped_crate_toml(manifest, "0.3.0", "0.3.1")

    def test_rejects_an_empty_unreleased_section(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            changelog = Path(temporary) / "CHANGELOG.md"
            with self.assertRaisesRegex(SystemExit, "nothing to release"):
                prepare_crate.require_nonempty_unreleased(changelog, "\nno entries here\n")

    def test_rolls_unreleased_into_a_dated_heading(self) -> None:
        prefix = "# Changelog\n\n"
        body = "\n### Added\n\n- one\n\n"
        suffix = ""
        rolled = prepare_crate.rolled_changelog(prefix, body, suffix, "0.3.1", "2026-09-22")
        self.assertIn(
            "## [Unreleased]\n\n## [0.3.1] - 2026-09-22\n### Added\n\n- one\n\n", rolled
        )

    def test_changelog_unreleased_body_splits_on_the_next_heading(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            changelog = Path(temporary) / "CHANGELOG.md"
            changelog.write_text(
                "# Changelog\n\n## [Unreleased]\n\n- pending\n\n## [0.3.0] - 2026-09-01\n\n- shipped\n",
                encoding="utf-8",
            )
            prefix, body, suffix = prepare_crate.changelog_unreleased_body(changelog)
            self.assertIn("- pending", body)
            self.assertNotIn("- shipped", body)
            self.assertIn("## [0.3.0]", suffix)


if __name__ == "__main__":
    unittest.main()
