"""Shared runner for the `scripts/probe_*.py` mutation probes (#282).

A probe declares its mutants as data and hands them to `main`, which owns
everything that must not differ between probes:

- Anchors fail loudly. Before a single test runs, every selected mutant's
  edits are applied in memory; an anchor that is missing or occurs more than
  once, or an edit that changes nothing, stops the probe with exit code 2 and
  names the mutant, the file and the anchor. A stale anchor never skips a
  mutant, so a refactor cannot quietly shrink what a probe verifies.
- Each test command runs in its own process group with a per-mutant timeout.
  A timeout kills the whole group (`os.killpg`) and counts as a kill, shown as
  `killed (timeout)`. An interrupt kills the group too.
- Mutated files are restored whatever happens (try/finally, also on
  KeyboardInterrupt and SIGTERM), and the restore is verified byte for byte.
- The summary is `N/M killed`; the exit code is non-zero when a mutant
  survived, a mutant did not compile where the probe forbids that, or an
  anchor failed.

Command line, shared by every probe:

    --list            print the mutant names and exit
    --only NAME       run only this mutant (repeatable; an exact name, or a
                      substring when no name matches exactly)
    --check-anchors   verify every anchor applies uniquely, run no tests
    --timeout SECS    override the per-mutant timeout

`scripts/check-probe-anchors.sh` runs `--check-anchors` on every probe as a
gate step.

Mutant tuples, in the shapes the probes already use (`path` and `tests` are
the probe-wide defaults given to `main`):

    (name, old, new)                          edit `path`
    (name, file, old, new)                    edit `file`, run `tests`
    (name, file, [(old, new), ...], target)   several edits in one file
    (name, file, old, new, target)            per-mutant target

A target is either one argument list (`["-p", "crate", "--test", "x"]`) or a
sequence of them; each runs as `command + target` until one fails. A `file`
is relative to `root` (default: the repository root) unless absolute.
"""
from __future__ import annotations

import argparse
import os
import signal
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

ROOT = Path(__file__).resolve().parents[1]
CARGO_TEST = ("cargo", "test", "-q")
EXIT_SURVIVED = 1
EXIT_ANCHOR = 2


@dataclass
class Mutant:
    name: str
    # (absolute path, old, new), applied in order.
    edits: list[tuple[Path, str, str]]
    # Full argv lists; the mutant is killed when any of them fails.
    commands: list[list[str]]


class AnchorError(Exception):
    pass


def _is_target(value) -> bool:
    """One argument list, as opposed to a sequence of them."""
    return isinstance(value, (list, tuple)) and (not value or isinstance(value[0], str))


def _targets(value) -> list[list[str]]:
    if value is None:
        return []
    if _is_target(value):
        return [list(value)]
    return [list(t) for t in value]


def normalize(raw, *, root: Path, path, tests, command: Sequence[str]) -> Mutant:
    """Read one mutant tuple in any of the shapes the module doc lists."""
    if isinstance(raw, Mutant):
        return raw
    name = raw[0]
    rest = list(raw[1:])
    if len(rest) == 2:  # (name, old, new) against the probe-wide path
        if path is None:
            raise TypeError(f"mutant {name!r} names no file and the probe gives no path")
        file, edits, target = path, [(rest[0], rest[1])], tests
    elif len(rest) == 3 and isinstance(rest[1], str):  # (name, file, old, new)
        file, edits, target = rest[0], [(rest[1], rest[2])], tests
    elif len(rest) == 3:  # (name, file, [(old, new), ...], target)
        file, edits, target = rest[0], list(rest[1]), rest[2]
    elif len(rest) == 4:  # (name, file, old, new, target)
        file, edits, target = rest[0], [(rest[1], rest[2])], rest[3]
    else:
        raise TypeError(f"mutant {name!r} has an unknown shape ({len(raw)} fields)")
    file = Path(file)
    if not file.is_absolute():
        file = root / file
    commands = [[*command, *t] for t in _targets(target)]
    if not commands:
        raise TypeError(f"mutant {name!r} has no test command")
    return Mutant(name, [(file, old, new) for old, new in edits], commands)


def _short(text: str, limit: int = 160) -> str:
    text = repr(text)
    return text if len(text) <= limit else text[: limit - 3] + "..."


def mutated_texts(mutant: Mutant, originals: dict[Path, str]) -> dict[Path, str]:
    """The mutated content of each file the mutant touches, or AnchorError."""
    texts: dict[Path, str] = {}
    for path, old, new in mutant.edits:
        if path not in texts:
            if path not in originals:
                if not path.is_file():
                    raise AnchorError(f"mutant {mutant.name!r}: file {path} does not exist")
                originals[path] = path.read_text()
            texts[path] = originals[path]
        count = texts[path].count(old)
        if count != 1:
            what = "not found" if count == 0 else f"found {count} times, not unique"
            raise AnchorError(
                f"mutant {mutant.name!r}: anchor {what} in {_rel(path)}: {_short(old)}"
            )
        if old == new:
            raise AnchorError(f"mutant {mutant.name!r}: edit in {_rel(path)} changes nothing")
        texts[path] = texts[path].replace(old, new, 1)
    return texts


def _rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def run_command(argv: list[str], *, cwd: Path, timeout: float, env) -> tuple[int | None, str]:
    """Exit code (None on timeout) and the combined output of one command.

    The command runs in its own session, so a timeout or an interrupt can
    kill everything it started (cargo, rustc, the test binaries) at once.
    """
    with tempfile.TemporaryFile(mode="w+", errors="replace") as out:
        process = subprocess.Popen(
            argv, cwd=cwd, env=env, stdout=out, stderr=subprocess.STDOUT,
            stdin=subprocess.DEVNULL, start_new_session=True,
        )
        try:
            code = process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            _kill_group(process)
            code = None
        except BaseException:
            _kill_group(process)
            raise
        out.seek(0)
        return code, out.read()


def _kill_group(process: subprocess.Popen) -> None:
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    process.wait()


def _compile_failed(output: str) -> bool:
    return "could not compile" in output or "error[E" in output


def _select(mutants: list[Mutant], only: list[str]) -> list[Mutant]:
    if not only:
        return mutants
    chosen: list[Mutant] = []
    for want in only:
        exact = [m for m in mutants if m.name == want]
        hits = exact or [m for m in mutants if want in m.name]
        if not hits:
            sys.exit(f"--only {want!r} matches no mutant (see --list)")
        chosen += [m for m in hits if m not in chosen]
    return [m for m in mutants if m in chosen]


def _raise_interrupt(signum, frame):
    raise KeyboardInterrupt(f"signal {signum}")


def main(
    mutants: Iterable,
    *,
    tests=None,
    path=None,
    root: Path | str | None = None,
    command: Sequence[str] = CARGO_TEST,
    cwd: Path | str | None = None,
    timeout: float = 1200,
    env: dict[str, str] | None = None,
    baseline: bool = False,
    reject_uncompiled: bool = False,
    argv: Sequence[str] | None = None,
) -> None:
    """Run the probe and exit.

    `tests` is the default target for mutants that name none; `command` is
    prefixed to every target (pass `command=()` for full argv lists);
    `baseline=True` runs the commands on the unmutated tree first and after
    the restore, and refuses to probe a red tree; `reject_uncompiled=True`
    counts a mutant that does not compile as BROKEN (a failure) rather than
    as killed.
    """
    parser = argparse.ArgumentParser(description=(sys.modules["__main__"].__doc__ or "").split("\n")[0])
    parser.add_argument("--list", action="store_true", help="print the mutant names and exit")
    parser.add_argument("--only", action="append", default=[], metavar="NAME",
                        help="run only this mutant (repeatable; exact name or substring)")
    parser.add_argument("--check-anchors", action="store_true",
                        help="verify every anchor applies uniquely; run no tests")
    parser.add_argument("--timeout", type=float, default=None, help="per-mutant timeout in seconds")
    args = parser.parse_args(argv)

    root = Path(root) if root is not None else ROOT
    cwd = Path(cwd) if cwd is not None else ROOT
    timeout = args.timeout if args.timeout is not None else timeout
    every = [normalize(m, root=root, path=path, tests=tests, command=command) for m in mutants]
    names = [m.name for m in every]
    duplicates = sorted({n for n in names if names.count(n) > 1})
    if duplicates:
        print(f"ANCHOR  duplicate mutant names: {duplicates}", file=sys.stderr)
        sys.exit(EXIT_ANCHOR)

    if args.list:
        for m in every:
            print(m.name)
        sys.exit(0)

    selected = _select(every, args.only)

    # Every anchor, before any test runs: a stale one stops the probe.
    originals: dict[Path, str] = {}
    failures = []
    for m in selected:
        try:
            mutated_texts(m, originals)
        except AnchorError as error:
            failures.append(str(error))
    if failures:
        for line in failures:
            print(f"ANCHOR  {line}", file=sys.stderr)
        print(f"{len(failures)}/{len(selected)} mutants have a stale anchor", file=sys.stderr)
        sys.exit(EXIT_ANCHOR)
    if args.check_anchors:
        print(f"anchors ok: {len(selected)} mutants")
        sys.exit(0)

    signal.signal(signal.SIGTERM, _raise_interrupt)
    survivors: list[str] = []
    broken: list[str] = []

    def run_all(mutant: Mutant) -> tuple[str, str, list[str] | None]:
        """Status of one run, stopping at the first red command, with that
        command's output and argv."""
        for argv_ in mutant.commands:
            code, output = run_command(argv_, cwd=cwd, timeout=timeout, env=env)
            if code is None:
                return "timeout", output, argv_
            if code != 0:
                return ("uncompiled" if _compile_failed(output) else "failed"), output, argv_
        return "passed", "", None

    def check_baseline(when: str) -> None:
        commands = {tuple(c) for m in selected for c in m.commands}
        for argv_ in sorted(commands):
            code, output = run_command(list(argv_), cwd=cwd, timeout=timeout, env=env)
            if code != 0:
                print(output[-3000:])
                print(f"{when}: FAIL ({' '.join(argv_)}) -- fix the suite before probing")
                sys.exit(EXIT_SURVIVED)
        print(f"{when}: PASS", flush=True)

    try:
        if baseline:
            check_baseline("baseline")
        for current in selected:
            texts = mutated_texts(current, originals)
            started = time.monotonic()
            try:
                for file, text in texts.items():
                    file.write_text(text)
                status, detail, killer = run_all(current)
            finally:
                for file in texts:
                    file.write_text(originals[file])
            seconds = time.monotonic() - started
            if status == "timeout":
                label = "killed (timeout)"
            elif status == "uncompiled":
                label = "BROKEN (does not compile)" if reject_uncompiled else "killed (compile error)"
            elif status == "failed":
                label = "killed"
            else:
                label = "SURVIVED"
            by = f" by {' '.join(killer)}" if killer and len(current.commands) > 1 else ""
            print(f"{label:26} {current.name}  [{seconds:.0f}s{by}]", flush=True)
            if status == "passed":
                survivors.append(current.name)
            elif status == "uncompiled" and reject_uncompiled:
                broken.append(current.name)
                print("\n".join(detail.splitlines()[-20:]), flush=True)
    finally:
        for file, text in originals.items():
            if file.read_text() != text:
                file.write_text(text)
        unrestored = [_rel(f) for f, t in originals.items() if f.read_text() != t]
        if unrestored:
            print(f"RESTORE FAILED: {unrestored}", file=sys.stderr)
            sys.exit(EXIT_ANCHOR)

    if baseline:
        check_baseline("restored")
    killed = len(selected) - len(survivors) - len(broken)
    print(f"{killed}/{len(selected)} killed")
    if survivors:
        print("SURVIVED: " + ", ".join(survivors))
    if broken:
        print("BROKEN (a mutant that does not compile proves nothing): " + ", ".join(broken))
    sys.exit(EXIT_SURVIVED if survivors or broken else 0)
