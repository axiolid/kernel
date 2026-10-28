#!/usr/bin/env bash
# Mutation-verify `cargo xtask architecture check`.
#
# A gate that has never failed is decoration with a green light. Each mutation
# below is a real architectural violation of the kind a hurried manifest edit
# would introduce. The gate must go RED for every one of them, and stay GREEN
# for the comment-only decoy.
#
# The script REFUSES to report a result when a mutation did not actually land
# (diff against the backup), because an unapplied patch and a blind gate look
# identical from the outside.
set -uo pipefail
export PATH="$HOME/.cargo/bin:$PATH"
cd "$(dirname "$0")/.."

GATE=(cargo xtask architecture check)
BAK="${TMPDIR:-/tmp}/layering_mut.$$.bak"
fail=0

# Content hash of every source/manifest in `crates/`, so we can prove the probe
# reverted exactly what it changed regardless of unrelated work in the tree.
snapshot() {
    find crates -type f \( -name '*.rs' -o -name '*.toml' \) \
        -exec md5sum {} + 2>/dev/null | sort -k2
}
BEFORE="$(snapshot)"

# Run the gate; echo "GREEN" or "RED".
run_gate() {
    if "${GATE[@]}" >"${TMPDIR:-/tmp}/layering_mut_out.$$.txt" 2>&1; then echo GREEN; else echo RED; fi
}

# mutate <label> <manifest> <line-to-append-to-[dependencies]> <expected>
mutate() {
    local label="$1" manifest="$2" line="$3" expect="$4"
    cp "$manifest" "$BAK"
    # Insert directly under the [dependencies] header.
    python3 - "$manifest" "$line" <<'PY'
import sys
path, line = sys.argv[1], sys.argv[2]
src = open(path).read()
assert "[dependencies]" in src, f"no [dependencies] table in {path}"
open(path, "w").write(src.replace("[dependencies]", "[dependencies]\n" + line, 1))
PY
    if diff -q "$manifest" "$BAK" >/dev/null; then
        echo "  $label: MUTATION DID NOT APPLY -- result would be meaningless"
        cp "$BAK" "$manifest"
        fail=1
        return
    fi
    local got; got=$(run_gate)
    cp "$BAK" "$manifest"
    if [ "$got" = "$expect" ]; then
        printf '  %-58s %s (expected %s)  ok\n' "$label" "$got" "$expect"
    else
        printf '  %-58s %s (expected %s)  MISS\n' "$label" "$got" "$expect"
        fail=1
    fi
}

# mutate_source <label> <rust-source> <line-to-append> <expected>
mutate_source() {
    local label="$1" source="$2" line="$3" expect="$4"
    cp "$source" "$BAK"
    printf '\n%s\n' "$line" >> "$source"
    if diff -q "$source" "$BAK" >/dev/null; then
        echo "  $label: MUTATION DID NOT APPLY -- result would be meaningless"
        cp "$BAK" "$source"
        fail=1
        return
    fi
    local got; got=$(run_gate)
    cp "$BAK" "$source"
    if [ "$got" = "$expect" ]; then
        printf '  %-58s %s (expected %s)  ok\n' "$label" "$got" "$expect"
    else
        printf '  %-58s %s (expected %s)  MISS\n' "$label" "$got" "$expect"
        fail=1
    fi
}

echo "=== baseline ==="
base=$(run_gate)
printf '  %-58s %s\n' "unmutated tree" "$base"
[ "$base" = GREEN ] || { echo "baseline is not green; fix that before mutating"; exit 1; }

echo "=== mutations ==="
G=crates

# 1. The seam reversed: geometry reaching back into IFC.
mutate "axiolid-mesh depends on ifc-model" \
    "$G/representations/discrete/mesh/Cargo.toml" "ifc-model.workspace = true" RED

# 2. Tier inversion: a representation crate pulling in an algorithm crate.
mutate "axiolid-mesh representation depends on dispatch execution" \
    "$G/representations/discrete/mesh/Cargo.toml" 'axiolid-dispatch.workspace = true' RED

# 3. The dependency is dev-only, but it is still absent from the package's
#    exact allowlist. Dev edges may point upward for integration tests only when
#    they are explicitly declared in architecture metadata.
cp "$G/foundation/core/Cargo.toml" "$BAK"
printf '\n[dev-dependencies]\naxiolid-mesh.workspace = true\n' >> "$G/foundation/core/Cargo.toml"
if diff -q "$G/foundation/core/Cargo.toml" "$BAK" >/dev/null; then
    echo "  axiolid-core dev-depends on axiolid-mesh: MUTATION DID NOT APPLY"; fail=1
else
    got=$(run_gate); cp "$BAK" "$G/foundation/core/Cargo.toml"
    if [ "$got" = RED ]; then
        printf '  %-58s %s (expected RED)  ok\n' "axiolid-core dev-depends on axiolid-mesh" "$got"
    else
        printf '  %-58s %s (expected RED)  MISS\n' "axiolid-core dev-depends on axiolid-mesh" "$got"; fail=1
    fi
fi

# 4. A new crate appearing without workspace registration or architecture metadata.
mkdir -p "$G/unregistered/src"
cat > "$G/unregistered/Cargo.toml" <<'TOML'
[package]
name = "axiolid-untiered"
version.workspace = true
edition.workspace = true

[dependencies]
axiolid-core.workspace = true
TOML
echo "// mutation probe" > "$G/unregistered/src/lib.rs"
got=$(run_gate)
rm -rf "$G/unregistered"
if [ "$got" = RED ]; then
    printf '  %-58s %s (expected RED)  ok\n' "new unregistered package manifest" "$got"
else
    printf '  %-58s %s (expected RED)  MISS\n' "new unregistered package manifest" "$got"; fail=1
fi

# 5. Protobuf transport vocabulary leaking into executable production source.
mutate_source "axiolid-mesh code names protobuf transport" \
    "$G/representations/discrete/mesh/src/lib.rs" 'const PROTOBUF_TRANSPORT: &str = "protobuf";' RED

# 6. Rust Protobuf implementation leaking into neutral production source.
mutate_source "axiolid-mesh source imports prost" \
    "$G/representations/discrete/mesh/src/lib.rs" "use prost::Message;" RED

# 7. Renaming a dependency key must not hide its package identity.
mutate "axiolid-mesh aliases prost dependency as p" \
    "$G/representations/discrete/mesh/Cargo.toml" \
    'p = { package = "prost", version = "0.13" }' RED

# 8. Comments are explanatory prose, not executable coupling.
mutate_source "Rust comment names protobuf/prost (must NOT trip)" \
    "$G/representations/discrete/mesh/src/lib.rs" \
    "// protobuf and prost documentation-only decoy" GREEN
mutate "COMMENTED-OUT ifc-model dep (must NOT trip)" \
    "$G/representations/discrete/mesh/Cargo.toml" "# ifc-model.workspace = true" GREEN

echo "=== restored ==="
printf '  %-58s %s\n' "tree after restore" "$(run_gate)"

# Verify THIS SCRIPT reverted its own mutations -- not that the tree is
# pristine. Demanding a clean worktree conflates "the probe leaked a mutation"
# with "the developer has unrelated work in progress", and the second is the
# normal case: you run this right after editing the gate or adding a crate.
# Compare against the snapshot taken at startup instead.
if [ "$(snapshot)" = "$BEFORE" ]; then
    echo "  crates identical to pre-probe state"
else
    echo "  DIRTY -- restore failed; probe mutations leaked:"
    diff <(printf '%s\n' "$BEFORE") <(snapshot) | sed 's/^/    /'
    fail=1
fi

echo
[ "$fail" -eq 0 ] && echo "MUTATION MATRIX PASSED" || echo "MUTATION MATRIX FAILED"
exit "$fail"
