#!/usr/bin/env bash
# Mutation probe: does the SONAME check actually fail without a SONAME?
# A guard that cannot fail is decoration (kernel#70).
set -uo pipefail
cd "$(dirname "$0")/.."
BUILD_RS=crates/facade/axiolid-capi/build.rs
BACKUP="${TMPDIR:-/tmp}/build.rs.bak"
target_dir="${CARGO_TARGET_DIR:-target}"
fail=0

check() {
  printf '  %-56s' "$1"
  if [ "$2" = "$3" ]; then echo "ok"; else echo "MISS (got=$2 want=$3)"; fail=1; fi
}

cp "$BUILD_RS" "$BACKUP" || exit 1
restore() { cp "$BACKUP" "$BUILD_RS"; rm -f "$BACKUP"; }
trap restore EXIT

echo "=== baseline: SONAME present, check passes ==="
cargo +1.88.0 build --release -p axiolid-capi >/dev/null 2>&1
have="$(readelf -d "${target_dir}/release/libaxiolid_capi.so" | grep -c SONAME)"
check "cdylib carries a SONAME" "$have" "1"

echo "=== mutation: drop the soname link arg ==="
cat > "$BUILD_RS" <<'RS'
fn main() {
    println!("cargo:rerun-if-changed=build.rs");
}
RS
touch crates/facade/axiolid-capi/src/lib.rs
cargo +1.88.0 build --release -p axiolid-capi >/dev/null 2>&1
have="$(readelf -d "${target_dir}/release/libaxiolid_capi.so" | grep -c SONAME)"
check "mutated artifact really loses its SONAME" "$have" "0"

if bash scripts/check-capi.sh >/dev/null 2>&1; then
  check "check-capi.sh rejects a soname-less cdylib" "accepted" "rejected"
else
  check "check-capi.sh rejects a soname-less cdylib" "rejected" "rejected"
fi

restore
trap - EXIT
touch crates/facade/axiolid-capi/src/lib.rs
cargo +1.88.0 build --release -p axiolid-capi >/dev/null 2>&1
have="$(readelf -d "${target_dir}/release/libaxiolid_capi.so" | grep -c SONAME)"
echo "=== restored ==="
check "SONAME back after restore" "$have" "1"

echo
[ "$fail" -eq 0 ] && echo "SONAME MATRIX PASSED" || echo "SONAME MATRIX FAILED"
exit "$fail"
