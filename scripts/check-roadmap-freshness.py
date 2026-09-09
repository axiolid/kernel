#!/usr/bin/env python3
"""Keep docs/ROADMAP.md from drifting back into a status page.

The roadmap explains *ordering and reasoning*. The moment it starts carrying
per-item status, it starts going stale — the exact failure this gate exists to
prevent.

Checks:
  1. No task checkboxes. A checklist here duplicates the project board.
  2. No progress words ("now:", "next:", "in progress", "completed") used as
     section headings, which imply a status this page cannot keep current.
  3. Every milestone that exists on GitHub is mentioned, so a new milestone
     cannot be invisible here. (Skipped without network/gh.)
  4. The pointer block to the board/milestones/capabilities is intact.

Run: python3 scripts/check-roadmap-freshness.py
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

ROADMAP = Path(__file__).resolve().parent.parent / "docs" / "ROADMAP.md"

# Headings that assert a state this page cannot keep current.
STALE_HEADINGS = re.compile(
    r"^#{2,3}\s+.*\b(now|next|then|in progress|completed|current status)\b",
    re.IGNORECASE | re.MULTILINE,
)
LIST_OR_REF = re.compile(r"(?m)^[ \t]*[-*+][ \t]|#\d+|/issues/")
CHECKBOX = re.compile(r"^\s*[-*]\s+\[[ xX]\]", re.MULTILINE)
REQUIRED_POINTERS = (
    "orgs/axiolid/projects/1",
    "axiolid/kernel/milestones",
    "capabilities",
)


def milestones_on_github() -> list[str] | None:
    """Milestone objects, or None when GitHub is unreachable."""
    try:
        out = subprocess.run(
            ["gh", "api", "repos/axiolid/kernel/milestones",
             "-X", "GET", "-f", "state=all", "--paginate"],
            capture_output=True, text=True, timeout=30,
        )
        if out.returncode != 0:
            return None
        return json.loads(out.stdout)
    except Exception:
        return None


def milestone_description_problems() -> list[str]:
    """A milestone description states one broad goal, nothing more.

    Checkboxes there render disabled and can never be ticked, and any
    restatement of the criteria duplicates the milestone's own issue list,
    which GitHub already tracks with a progress bar.
    """
    try:
        out = subprocess.run(
            ["gh", "api", "repos/axiolid/kernel/milestones",
             "-X", "GET", "-f", "state=all", "--paginate"],
            capture_output=True, text=True, timeout=30, check=True,
        ).stdout
    except Exception:
        return []
    bad = []
    for m in json.loads(out):
        d = m.get("description") or ""
        # A milestone states one broad goal; the issue list IS the
        # criteria and GitHub tracks its progress. Any list here
        # duplicates that and goes stale.
        if LIST_OR_REF.search(d):
            bad.append(m["title"])
    return bad


def main() -> int:
    if not ROADMAP.exists():
        print(f"roadmap: {ROADMAP} not found")
        return 1

    text = ROADMAP.read_text(encoding="utf-8")
    problems: list[str] = []

    for match in CHECKBOX.finditer(text):
        line = text[: match.start()].count("\n") + 1
        problems.append(
            f"line {line}: task checkbox — per-item status belongs on the "
            f"project board, not here"
        )

    for match in STALE_HEADINGS.finditer(text):
        line = text[: match.start()].count("\n") + 1
        heading = match.group(0).strip()
        problems.append(
            f"line {line}: heading {heading!r} asserts progress state this "
            f"page cannot keep current"
        )

    for pointer in REQUIRED_POINTERS:
        if pointer not in text:
            problems.append(
                f"missing pointer to {pointer!r} — readers must be sent to the "
                f"live source"
            )

    milestones = milestones_on_github()
    if milestones is None:
        print("roadmap: skipping milestone coverage (gh unavailable)")
    else:
        for milestone in milestones:
            # Require the milestone's own LINK, not its prose. Matching title
            # text anywhere in the page is vacuous: the words recur in the
            # surrounding argument, so deleting a whole row still passed.
            # Verified by mutation -- removing a row must fail this check.
            link = f"/milestone/{milestone['number']})"
            if link not in text:
                problems.append(
                    f"milestone {milestone['title']!r} exists on GitHub but "
                    f"has no row linking to {link} here"
                )

        # And the reverse: a row may not link to a milestone that does not
        # exist. Without this, repointing a row at /milestone/999 silently
        # passes -- the row looks present while linking nowhere real.
        live_numbers = {m["number"] for m in milestones}
        for linked in re.findall(r"/milestone/(\d+)\)", text):
            if int(linked) not in live_numbers:
                problems.append(
                    f"roadmap links to /milestone/{linked} which does not "
                    f"exist on GitHub"
                )

    for title in milestone_description_problems():
        problems.append(
            f"milestone {title!r} carries a list or issue reference in its "
            f"description. The milestone's own issue list is the criteria and "
            f"GitHub tracks its progress; keep the description to one broad goal."
        )

    if problems:
        print(f"roadmap: {len(problems)} problem(s):")
        for problem in problems:
            print(f"- {problem}")
        return 1

    print("roadmap freshness: ok")
    return 0


if __name__ == "__main__":
    sys.exit(main())
