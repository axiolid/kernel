#!/usr/bin/env bash
# Deterministic instruction-count regression benchmarks.
#
# Skips cleanly when the tooling is absent instead of failing: a developer
# without valgrind should still be able to run the gate, and a benchmark that
# cannot run is not a regression. CI installs the tooling, so the gate is real
# there.
set -euo pipefail

cd "$(dirname "$0")/.."

if ! command -v valgrind > /dev/null 2>&1; then
    echo "SKIP: valgrind not installed (apt-get install valgrind)"
    exit 0
fi

declared="$(grep -m1 '^iai-callgrind' tools/benchmark/Cargo.toml | sed 's/.*"\(.*\)".*/\1/')"

if ! command -v iai-callgrind-runner > /dev/null 2>&1; then
    echo "SKIP: iai-callgrind-runner not installed"
    echo "      cargo install iai-callgrind-runner --version ${declared}"
    exit 0
fi

# The runner must match the library version or it refuses to run. It has no
# plain --version output (it reports a diagnostic instead), so the mismatch is
# detected from the benchmark run itself rather than by parsing a version
# string that does not exist.
exec cargo bench -p axiolid-benchmark --bench regression "$@"
