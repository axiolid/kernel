#!/usr/bin/env python3
"""Fail when the public API breaks without a version bump that admits it.

A policy nobody can check is not in force (kernel#34). This wraps
cargo-semver-checks, whose baseline is the crate actually published on
crates.io rather than a snapshot committed in-tree -- an in-tree
baseline drifts, and a drifted baseline reports success for a break.

Pre-1.0 Cargo treats the MINOR field as the major, so 0.2.x -> 0.3.0 is
the breaking step and 0.2.0 -> 0.2.1 is not. The check therefore only
has to answer one question: does the working tree break the published
API without the workspace version having already moved past it?
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def workspace_version() -> str:
    """The single version every publishable crate in this workspace carries."""
    meta = json.loads(
        subprocess.run(
            ["cargo", "metadata", "--format-version", "1", "--no-deps", "--locked"],
            cwd=ROOT,
            check=True,
            text=True,
            stdout=subprocess.PIPE,
        ).stdout
    )
    versions = {p["version"] for p in meta["packages"] if p.get("publish") != []}
    if len(versions) != 1:
        raise SystemExit(
            f"expected one workspace version, found {sorted(versions)}"
        )
    return versions.pop()


def publishable() -> list[str]:
    """Crates a consumer can depend on, so crates whose API we owe."""
    meta = json.loads(
        subprocess.run(
            ["cargo", "metadata", "--format-version", "1", "--no-deps", "--locked"],
            cwd=ROOT,
            check=True,
            text=True,
            stdout=subprocess.PIPE,
        ).stdout
    )
    return sorted(p["name"] for p in meta["packages"] if p.get("publish") != [])


def published_versions(name: str) -> list[str]:
    """Released versions from the sparse index, newest last.

    The sparse index is used rather than the crates.io web API because
    the latter rejects unfamiliar user agents with a 403, which would
    turn a network policy into a silent gate skip.
    """
    url = f"https://index.crates.io/{name[:2]}/{name[2:4]}/{name}"
    if len(name) < 4:
        url = f"https://index.crates.io/{len(name)}/{name}"
    request = urllib.request.Request(
        url, headers={"User-Agent": "axiolid-semver-gate/1.0"}
    )
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            body = response.read().decode("utf-8")
    except urllib.error.HTTPError as error:
        if error.code == 404:
            return []
        raise
    out = []
    for line in body.splitlines():
        if not line.strip():
            continue
        entry = json.loads(line)
        if not entry.get("yanked"):
            out.append(entry["vers"])
    return out


def parse(version: str) -> tuple[int, int, int]:
    major, minor, patch = (int(part) for part in version.split(".")[:3])
    return major, minor, patch


def bump_admits_breakage(current: str, published: str) -> bool:
    """Has the version already moved in a way that allows a break?

    Cargo's compatibility rule for 0.x is the MINOR field, so 0.2.1 is
    compatible with 0.2.0 and 0.3.0 is not. Post-1.0 it is MAJOR.
    """
    cur = parse(current)
    pub = parse(published)
    if pub[0] == 0:
        return (cur[0], cur[1]) > (pub[0], pub[1])
    return cur[0] > pub[0]


def newest(versions: list[str]) -> str:
    return max(versions, key=parse)


def run_semver_checks(crates, baseline):
    """Check the given crates against an explicitly pinned baseline.

    The baseline must be named. Left to itself cargo-semver-checks picks
    the newest release, which here equals the working-tree version, so it
    compares the tree against itself and can never report a break.
    """
    args = ["cargo", "semver-checks", "check-release", "--baseline-version", baseline]
    for name in crates:
        args += ["--package", name]
    process = subprocess.run(
        args, cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
    )
    return process.returncode, process.stdout


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--explain", action="store_true")
    args = parser.parse_args()

    if shutil.which("cargo-semver-checks") is None:
        # Fail loudly rather than skip. A gate that silently disappears when
        # a tool is missing is worse than no gate: it reports success.
        print("cargo-semver-checks is not installed; the breaking-change")
        print("gate cannot run. Install it with:")
        print("    cargo install cargo-semver-checks --version 0.44.0 --locked")
        return 1

    current = workspace_version()
    names = publishable()

    checkable: list[str] = []
    skipped: list[str] = []
    admitted: list[str] = []
    baselines: dict[str, list[str]] = {}
    for name in names:
        versions = published_versions(name)
        if not versions:
            skipped.append(name)
            continue
        # Compare against the newest release STRICTLY BELOW the working
        # version. The equal-version case is the tree against itself, which
        # cannot fail and would make the gate decorative.
        prior = [v for v in versions if parse(v) < parse(current)]
        if not prior:
            skipped.append(name)
            continue
        baseline = newest(prior)
        if bump_admits_breakage(current, baseline):
            admitted.append(f"{name} {baseline} -> {current}")
            continue
        baselines.setdefault(baseline, []).append(name)
        checkable.append(name)

    if args.explain:
        print(f"workspace version: {current}")
        print(f"checked against crates.io : {len(checkable)}")
        print(f"unpublished, no baseline  : {len(skipped)}")
        print(f"bump already admits break : {len(admitted)}")
        for entry in admitted:
            print(f"    {entry}")
        return 0

    if not checkable:
        # Every crate is either unpublished or already past its baseline.
        # Nothing can break compatibility that nobody can depend on yet.
        print("semver: nothing to check (no published baseline in scope)")
        return 0

    # Crates can sit on different baselines, so group by baseline and run
    # one pass per group rather than assuming a single shared version.
    failed = False
    collected = []
    for baseline, group in sorted(baselines.items()):
        code, output = run_semver_checks(group, baseline)
        collected.append(output)
        if code != 0:
            failed = True
    if failed:
        print("\n".join(collected))
        print("")
        print("A breaking change needs a version bump that admits it.")
        print("See docs/contributing/breaking-changes.md")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
