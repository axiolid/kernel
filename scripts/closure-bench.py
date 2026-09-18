"""Measure what each closure profile costs to build (kernel#12).

ADR 0036 recorded a spot check: cold wall time and target/ size for
two profiles, best of three. This makes it repeatable and widens it
to every declared profile.

Each profile gets its OWN CARGO_TARGET_DIR so no profile is measured
against another one warm cache. The noise floor is measured first,
by rebuilding one profile repeatedly, so a reader can tell a real
separation from run-to-run jitter.
"""

import argparse
import json
import os
import shutil
import statistics
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CONSUMERS = ROOT / "tests" / "consumers"
BENCH_ROOT = Path("/mnt/backup/build-cache/closure-bench")


def profiles():
    """Every declared closure profile, discovered not hardcoded.

    A profile added later is measured automatically; a list here would
    silently omit it and the table would look complete anyway.
    """
    found = []
    for path in sorted(CONSUMERS.iterdir()):
        if (path / "Cargo.toml").is_file():
            found.append(path.name)
    return found


def env_for(target_dir):
    """A clean environment: the shared shell leaks cargo vars.

    RUSTUP_TOOLCHAIN is dropped so rust-toolchain.toml decides, and
    every arm is built by the same pinned compiler. A leaked value
    would silently compare two toolchains.
    """
    env = dict(os.environ)
    for key in ("CARGO_TARGET_DIR", "CARGO_HOME", "GIT_DIR",
                "GIT_WORK_TREE", "RUSTUP_TOOLCHAIN", "CARGO_NET_OFFLINE"):
        env.pop(key, None)
    env["CARGO_TARGET_DIR"] = str(target_dir)
    return env


def dir_size_mb(path):
    total = 0
    for entry in path.rglob("*"):
        try:
            if entry.is_file() and not entry.is_symlink():
                total += entry.stat().st_size
        except OSError:
            continue
    return total / (1024 * 1024)


def package_count(manifest, env):
    """Resolved packages in this profile's closure.

    The number ADR 0036 argues about, read from cargo rather than
    counted by hand.
    """
    out = subprocess.run(
        ["cargo", "metadata", "--format-version", "1",
         "--manifest-path", str(manifest)],
        env=env, cwd=ROOT, stdout=subprocess.PIPE, text=True, check=True)
    return len(json.loads(out.stdout)["packages"])


def cold_build(name, release, env_target):
    """One COLD build: the target dir is wiped first.

    Without the wipe this measures an incremental no-op and every
    profile looks identically fast.
    """
    manifest = CONSUMERS / name / "Cargo.toml"
    target = env_target / name
    if target.exists():
        shutil.rmtree(target)
    target.mkdir(parents=True, exist_ok=True)
    env = env_for(target)
    cmd = ["cargo", "build", "--quiet", "--manifest-path", str(manifest)]
    if release:
        cmd.append("--release")
    start = time.perf_counter()
    proc = subprocess.run(cmd, env=env, cwd=ROOT,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, text=True)
    elapsed = time.perf_counter() - start
    if proc.returncode != 0:
        return None, proc.stdout
    return elapsed, target


def noise_floor(name, reps, release, bench_root):
    """Rebuild ONE profile repeatedly; the spread is the floor.

    Any difference between two profiles smaller than this is not a
    result. Measured rather than assumed: this box is shared.
    """
    samples = []
    for _ in range(reps):
        elapsed, _ = cold_build(name, release, bench_root)
        if elapsed is None:
            return None
        samples.append(elapsed)
    median = statistics.median(samples)
    spread = (max(samples) - min(samples)) / median * 100.0
    return {"samples": samples, "median": median, "spread_pct": spread}


def measure(name, reps, release, bench_root):
    manifest = CONSUMERS / name / "Cargo.toml"
    env = env_for(bench_root / name)
    samples = []
    size = 0.0
    for _ in range(reps):
        elapsed, target = cold_build(name, release, bench_root)
        if elapsed is None:
            return {"profile": name, "error": str(target)[:200]}
        samples.append(elapsed)
        size = dir_size_mb(target)
    return {
        "profile": name,
        "packages": package_count(manifest, env),
        "median_s": statistics.median(samples),
        "min_s": min(samples),
        "max_s": max(samples),
        "target_mb": size,
        "samples": samples,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reps", type=int, default=3)
    parser.add_argument("--release", action="store_true")
    parser.add_argument("--only", default=None)
    parser.add_argument("--json", default=None)
    parser.add_argument("--bench-root", default=str(BENCH_ROOT))
    args = parser.parse_args()

    bench_root = Path(args.bench_root)
    bench_root.mkdir(parents=True, exist_ok=True)
    names = profiles()
    if args.only:
        wanted = set(args.only.split(","))
        names = [n for n in names if n in wanted]
    if not names:
        print("no closure profiles found")
        return 1

    mode = "release" if args.release else "dev"
    print(f"closure profiles: {len(names)}  reps={args.reps}  mode={mode}")

    floor = noise_floor(names[0], max(2, args.reps), args.release, bench_root)
    if floor is None:
        print(f"noise floor could not be measured: {names[0]} failed to build")
        return 1
    print(f"noise floor on {names[0]}: {floor['spread_pct']:.1f}% "
          f"of a {floor['median']:.2f}s median")
    print("differences below that are not results")
    print()

    rows = [measure(n, args.reps, args.release, bench_root) for n in names]

    print(f"{'profile':<30}{'pkgs':>6}{'median s':>11}{'target MB':>12}")
    for row in rows:
        if "error" in row:
            print(f"{row['profile']:<30}{'BUILD FAILED':>29}")
            continue
        print(f"{row['profile']:<30}{row['packages']:>6}"
              f"{row['median_s']:>11.2f}{row['target_mb']:>12.1f}")

    if args.json:
        payload = {"mode": mode, "reps": args.reps,
                   "noise_floor": floor, "profiles": rows}
        Path(args.json).write_text(json.dumps(payload, indent=2))
        print(f"\nraw samples written to {args.json}")
    return 0 if all("error" not in r for r in rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
