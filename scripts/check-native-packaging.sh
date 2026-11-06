#!/usr/bin/env bash
set -euo pipefail

python3 -m py_compile scripts/package-native.py scripts/verify-native-package.py scripts/test-native-cmake.py scripts/test-native-fetch.py scripts/verify-native-release-set.py
python3 -m unittest tests/native/test_native_packaging.py
python3 -m unittest tests/native/test_package_budget.py
if grep -R -E 'target-cpu=native|CMAKE_(C|CXX)_FLAGS' native; then
  printf '%s\n' 'native integration must not set host-specific codegen or consumer-global flags' >&2
  exit 1
fi

log="$(mktemp)"
trap 'rm -f "$log"' EXIT
if cmake -P tests/native/reject-mutable-ref.cmake >"$log" 2>&1; then
  printf 'mutable source ref refusal unexpectedly passed\n' >&2
  exit 1
fi
if ! python3 - "$log" <<'PY'
from pathlib import Path
import sys
if "immutable 40-hex GIT_COMMIT" not in Path(sys.argv[1]).read_text():
    raise SystemExit(1)
PY
then
  printf 'mutable source ref failed for the wrong reason\n' >&2
  exit 1
fi

python3 scripts/test-native-cmake.py --build-type Release --linkage SHARED --mutations
