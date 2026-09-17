#!/usr/bin/env python3
"""Every capability row must name evidence that actually exists.

The v1.0 milestone requires that no published capability row is
aspirational (kernel#35). The risk is silent drift: a row stays green
after the code behind it moves or is renamed, and nothing notices.

A row's evidence cell names crates in backticks. This checks that each
named `axiolid-*` crate is a real workspace member, so a rename or a
removal breaks the gate instead of quietly leaving a false claim.

What this deliberately does NOT do is judge whether the prose is true.
That is not mechanically decidable. It enforces the weaker property
the page itself promises: the evidence pointer resolves.
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PAGE = ROOT / "docs" / "capabilities.md"

# A markdown table row: leading pipe, at least three cells.
ROW = re.compile(r"^\|(?P<cells>.+)\|\s*$")
CRATE = re.compile(r"`(axiolid[a-z0-9-]*)`")
ADR = re.compile(r"\(([^)]*?/adr/(\d{4})-[a-z0-9-]+)\)")


def workspace_crates() -> set[str]:
    meta = json.loads(
        subprocess.run(
            ["cargo", "metadata", "--format-version", "1", "--no-deps", "--locked"],
            cwd=ROOT,
            check=True,
            text=True,
            stdout=subprocess.PIPE,
        ).stdout
    )
    return {p["name"] for p in meta["packages"]}


def main() -> int:
    if not PAGE.exists():
        print(f"capabilities: {PAGE} not found")
        return 1

    crates = workspace_crates()
    text = PAGE.read_text(encoding="utf-8")
    problems: list[str] = []
    rows = 0
    cited = 0

    for number, line in enumerate(text.splitlines(), start=1):
        match = ROW.match(line)
        if not match:
            continue
        cells = [c.strip() for c in match.group("cells").split("|")]
        if len(cells) < 3:
            continue
        if cells[0] in {"Capability", "---"} or set(cells[0]) <= {"-"}:
            continue
        rows += 1
        # Evidence prose contains '|' (feature lists, alternations), so the
        # cell split over-fragments it. Everything after the capability name
        # is evidence; searching the whole remainder avoids depending on a
        # cell count the prose controls.
        evidence = "|".join(cells[1:])

        named = CRATE.findall(evidence)
        for crate in named:
            if crate not in crates:
                problems.append(
                    f"line {number}: evidence names `{crate}`, which is not a "
                    f"workspace crate — the claim points at nothing"
                )

        for link, adr_number in ADR.findall(evidence):
            matches = list((ROOT / "docs" / "adr").glob(f"{adr_number}-*.md"))
            if not matches:
                problems.append(
                    f"line {number}: evidence cites ADR {adr_number}, which "
                    f"does not exist"
                )

        if named or ADR.search(evidence):
            cited += 1
            continue

        # A row with no crate and no ADR is prose asserting a capability
        # with nothing a reader can check.
        problems.append(
            f"line {number}: row {cells[0]!r} cites no crate and no ADR"
        )

    if problems:
        print("capabilities: evidence does not resolve")
        for problem in problems:
            print(f"  {problem}")
        print("")
        print("Every row must point at something a reader can open.")
        return 1

    print(f"capabilities: {rows} rows, {cited} with resolvable evidence")
    return 0


if __name__ == "__main__":
    sys.exit(main())
