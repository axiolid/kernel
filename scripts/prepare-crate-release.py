#!/usr/bin/env python3
"""Bump ONE crate's version and roll its own CHANGELOG.md.

This is the per-crate counterpart to prepare-release.py, which bumps the
whole workspace at once for coordinated major releases. A single crate can
publish a patch or minor (pre-1.0: only a patch, since Cargo's caret rule
treats a 0.x minor as breaking) without touching any other crate's version,
because internal dependents already carry a caret requirement
(`version = "0.3"` in [workspace.dependencies]) that the new version keeps
satisfying.

Usage:
    python3 scripts/prepare-crate-release.py --crate axiolid-core --release 0.3.1 [--check]

--check (default) validates the bump and changelog shape without writing.
Pass --write to apply it. A bump that would break the existing internal
dependency requirement is refused: that is a workspace-wide event and belongs
to prepare-release.py instead.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CARGO_TOML = ROOT / "Cargo.toml"

SEMVER_RE = re.compile(r"^(0|[1-9]\d*)\.(0|[1-9]\d*)\.(0|[1-9]\d*)$")
UNRELEASED_HEADING = "## [Unreleased]"


def parse_semver(value: str) -> tuple[int, int, int]:
    match = SEMVER_RE.match(value)
    if match is None:
        raise SystemExit(f"not a valid semantic version: {value!r}")
    return tuple(int(part) for part in match.groups())  # type: ignore[return-value]


def require_forward_bump(current: str, requested: str) -> None:
    if parse_semver(requested) <= parse_semver(current):
        raise SystemExit(
            f"requested release {requested} must be strictly greater than "
            f"the current crate version {current}"
        )


def require_compatible_bump(current: str, requested: str) -> None:
    """Refuse a bump that an existing internal caret requirement would reject.

    Cargo's caret rule on 0.x treats the minor as the breaking slot:
    `^0.3.0` matches 0.3.* but not 0.4.0. At 1.x+ the major is the breaking
    slot as usual. A crate whose own version moves past what its existing
    dependents require is a workspace-wide event -- every internal
    `version = "..."` pin in [workspace.dependencies] would need editing, and
    that is exactly the multi-crate republish this mechanism exists to avoid
    doing by accident. Route it through prepare-release.py instead, which
    edits every pin in one pass.
    """
    current_major, current_minor, _ = parse_semver(current)
    new_major, new_minor, _ = parse_semver(requested)
    if current_major == 0:
        compatible = new_major == 0 and new_minor == current_minor
    else:
        compatible = new_major == current_major
    if not compatible:
        raise SystemExit(
            f"{requested} is a breaking bump from {current} under Cargo's "
            "caret rule (0.x treats the minor as breaking; 1.x+ treats the "
            "major as breaking). A breaking bump changes what every internal "
            "dependent requires, so it is a workspace-wide event -- use "
            "prepare-release.py to bump every crate and pin together."
        )


def metadata() -> dict:
    result = subprocess.run(
        ["cargo", "metadata", "--format-version", "1", "--no-deps", "--locked"],
        cwd=ROOT,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    )
    return json.loads(result.stdout)


def find_crate(name: str) -> tuple[Path, str]:
    """Return (manifest_path, current_version) for a publishable crate name."""
    data = metadata()
    members = set(data["workspace_members"])
    for package in data["packages"]:
        if package["id"] not in members or package.get("publish") == []:
            continue
        if package["name"] == name:
            return Path(package["manifest_path"]), package["version"]
    raise SystemExit(f"no publishable workspace crate named {name!r}")


CRATE_VERSION_RE = re.compile(r'(?m)^version = "(\d+\.\d+\.\d+)"$')


def bumped_crate_toml(manifest: Path, current: str, version: str) -> str:
    text = manifest.read_text(encoding="utf-8")
    replaced, count = CRATE_VERSION_RE.subn(f'version = "{version}"', text, count=1)
    if count != 1:
        raise SystemExit(f"{manifest}: expected exactly one version line, found {count}")
    if f'version = "{current}"' not in text:
        raise SystemExit(
            f"{manifest}: did not contain the expected current version {current!r}; "
            "metadata and the manifest disagree"
        )
    return replaced


def changelog_unreleased_body(changelog: Path) -> tuple[str, str, str]:
    text = changelog.read_text(encoding="utf-8")
    start = text.find(UNRELEASED_HEADING)
    if start == -1:
        raise SystemExit(f"{changelog}: missing {UNRELEASED_HEADING!r} heading")
    body_start = start + len(UNRELEASED_HEADING)
    next_heading = re.search(r"(?m)^## \[", text[body_start:])
    body_end = body_start + next_heading.start() if next_heading else len(text)
    return text[:start], text[body_start:body_end], text[body_end:]


def require_nonempty_unreleased(changelog: Path, body: str) -> None:
    if not any(line.strip().startswith("- ") for line in body.splitlines()):
        raise SystemExit(
            f"{changelog}: {UNRELEASED_HEADING} has no '- ' entries; "
            "nothing to release, or an entry was forgotten"
        )


def rolled_changelog(prefix: str, body: str, suffix: str, version: str, today: str) -> str:
    dated_heading = f"## [{version}] - {today}"
    return f"{prefix}{UNRELEASED_HEADING}\n\n{dated_heading}{body}{suffix}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--crate", required=True, help="publishable crate name, e.g. axiolid-core")
    parser.add_argument("--release", required=True, help="target version, e.g. 0.3.1")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check", action="store_true", default=True, help="validate only (default)")
    mode.add_argument("--write", action="store_true", help="apply the bump and changelog roll")
    args = parser.parse_args()

    manifest, current = find_crate(args.crate)
    require_forward_bump(current, args.release)
    require_compatible_bump(current, args.release)

    changelog = manifest.parent / "CHANGELOG.md"
    if not changelog.exists():
        raise SystemExit(f"{changelog}: does not exist; run scaffold-crate-changelogs.py first")
    prefix, body, suffix = changelog_unreleased_body(changelog)
    require_nonempty_unreleased(changelog, body)

    today = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    new_changelog = rolled_changelog(prefix, body, suffix, args.release, today)
    new_toml = bumped_crate_toml(manifest, current, args.release)

    entries = sum(1 for line in body.splitlines() if line.strip().startswith("- "))
    print(f"CRATE_RELEASE_CHECK crate={args.crate} current={current} requested={args.release} date={today}")
    print(f"CRATE_RELEASE_CHECK unreleased_entries={entries}")

    if not args.write:
        print("CRATE_RELEASE_PREPARE=CHECK_ONLY (pass --write to apply)")
        return 0

    manifest.write_text(new_toml, encoding="utf-8")
    changelog.write_text(new_changelog, encoding="utf-8")
    print(f"CRATE_RELEASE_PREPARE=WRITE crate={args.crate} version={args.release}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
