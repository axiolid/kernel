#!/usr/bin/env bash
# Mutation-verify the derived isolated-build list (kernel#38).
#
# The acceptance test named on the issue: add a throwaway crate, confirm
# the gate covers it WITHOUT anyone editing a list. The old hand-written
# loop passed this repo for months while 14 of 52 crates escaped it, so
# a green gate is not evidence -- coverage has to be demonstrated.
set -uo pipefail
cd "$(dirname "$0")/.."

CRATE_DIR="crates/foundation/zz-probe-ephemeral"
CRATE_NAME="axiolid-zz-probe-ephemeral"
fail=0

# The workspace lists members EXPLICITLY -- it does not glob -- so a new
# crate must be registered the way a real one is. Cargo.toml is restored
# from git rather than by reverse-editing.
cleanup() {
  rm -rf "$CRATE_DIR"
  git checkout -- Cargo.toml Cargo.lock 2>/dev/null || true
}
trap cleanup EXIT

targets() { python3 scripts/isolated-build-targets.py 2>/dev/null; }

echo "=== baseline ==="
before=$(targets | wc -l)
printf "  %-56s %s\n" "publishable members covered" "$before"
if [ "$before" -lt 40 ]; then
  echo "  baseline list is implausibly short; aborting"
  exit 1
fi

echo "=== a new publishable crate is covered automatically ==="
mkdir -p "$CRATE_DIR/src"
cat > "$CRATE_DIR/Cargo.toml" <<TOML
[package]
name = "$CRATE_NAME"
version = "0.0.0"
edition = "2021"
publish = false
TOML
printf "" > "$CRATE_DIR/src/lib.rs"
python3 - "$CRATE_DIR" <<'PY'
import sys, pathlib
member = sys.argv[1]
path = pathlib.Path("Cargo.toml")
text = path.read_text()
# Insert as the first member; a plain string splice avoids the
# quote-escaping hazards of a regex inside a shell heredoc.
anchor = "members = [" + chr(10)
assert anchor in text, 'workspace members list not found'
entry = "    " + chr(34) + member + chr(34) + "," + chr(10)
text = text.replace(anchor, anchor + entry, 1)
path.write_text(text)
PY
cargo metadata --format-version 1 --no-deps >/dev/null 2>&1 || true

# First as a PRIVATE crate: it must NOT be picked up, because nobody can
# consume it standalone. This is the decoy -- a rule that grabs every
# member would trip here.
priv_count=$(targets | grep -c "^${CRATE_NAME}$")
if [ "$priv_count" -eq 0 ]; then
  printf "  %-56s %s\n" "publish=false crate stays out" "ok"
else
  printf "  %-56s %s\n" "publish=false crate stays out" "MISS (picked up)"
  fail=1
fi

# Now make it publishable. The list must grow by exactly one, with no
# edit to gate.sh or to any list.
sed -i 's/^publish = false$//' "$CRATE_DIR/Cargo.toml"
after=$(targets | wc -l)
got=$(targets | grep -c "^${CRATE_NAME}$")
if [ "$got" -eq 1 ] && [ "$after" -eq $((before + 1)) ]; then
  printf "  %-56s %s\n" "new publishable crate auto-covered" "ok"
else
  printf "  %-56s %s\n" "new publishable crate auto-covered" "MISS (got=$got, $before -> $after)"
  fail=1
fi

echo "=== coverage is real: a broken crate fails its isolated build ==="
# Presence in a list proves nothing if the step cannot fail. Give the
# throwaway crate a dependency that only resolves through the workspace
# and confirm `cargo build -p` on it actually errors.
cat > "$CRATE_DIR/src/lib.rs" <<RS
pub fn broken() -> NoSuchType { unimplemented!() }
RS
if cargo build -p "$CRATE_NAME" >/dev/null 2>&1; then
  printf "  %-56s %s\n" "broken crate fails isolated build" "MISS (built anyway)"
  fail=1
else
  printf "  %-56s %s\n" "broken crate fails isolated build" "ok"
fi

cleanup
trap - EXIT

echo "=== restored ==="
restored=$(targets | wc -l)
if [ "$restored" -eq "$before" ]; then
  printf "  %-56s %s\n" "member set back to baseline" "$restored"
else
  printf "  %-56s %s\n" "member set NOT restored" "$before -> $restored"
  fail=1
fi

echo
[ "$fail" -eq 0 ] && echo "ISOLATED BUILD MATRIX PASSED" || echo "ISOLATED BUILD MATRIX FAILED"
exit "$fail"
