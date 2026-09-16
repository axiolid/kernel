#!/usr/bin/env bash
# Mutation-verify the crate-naming rule in the architecture gate.
#
# A naming rule nobody can violate is a comment. Each mutation below is a
# real violation of the kind a new crate would introduce; the gate must go
# RED for every one, and stay GREEN for the decoys that only LOOK wrong.
#
# Refuses to report when a mutation did not land: an unapplied patch and a
# blind gate are indistinguishable from the outside.
set -uo pipefail
export PATH="$HOME/.cargo/bin:$PATH"
cd "$(dirname "$0")/.."

GATE=(cargo xtask architecture check)
BAK="${TMPDIR:-/tmp}/naming_mut.bak"
fail=0

snapshot() {
    find crates tools -type f -name '*.toml' -exec md5sum {} + 2>/dev/null | sort -k2
}
BEFORE="$(snapshot)"

run_gate() {
    if "${GATE[@]}" >"${TMPDIR:-/tmp}/naming_mut_out.txt" 2>&1; then echo GREEN; else echo RED; fi
}

# swap <label> <manifest> <old> <new> <expected>
swap() {
    local label="$1" manifest="$2" old="$3" new="$4" expect="$5"
    cp "$manifest" "$BAK"
    python3 - "$manifest" "$old" "$new" <<'PY'
import sys
path, old, new = sys.argv[1], sys.argv[2], sys.argv[3]
src = open(path).read()
assert src.count(old) == 1, f'{old!r} not unique in {path}'
open(path, 'w').write(src.replace(old, new, 1))
PY
    if diff -q "$manifest" "$BAK" >/dev/null; then
        echo "  $label: MUTATION DID NOT APPLY -- result would be meaningless"
        cp "$BAK" "$manifest"; fail=1; return
    fi
    local got; got=$(run_gate)
    cp "$BAK" "$manifest"
    if [ "$got" = "$expect" ]; then
        printf '  %-58s %s (expected %s)  ok\n' "$label" "$got" "$expect"
    else
        printf '  %-58s %s (expected %s)  MISS\n' "$label" "$got" "$expect"; fail=1
    fi
}

# Every mutation below edits `[package.metadata.axiolid]`, never `name`.
# Renaming a crate breaks its dependents' manifests, so cargo fails to
# RESOLVE and the gate never runs -- a red light that proves nothing about
# this rule. Mutating role/domain leaves the graph intact, so a RED result
# is attributable to the naming check itself. Each case is asserted to
# mention the naming rule, not merely to fail.

echo "=== baseline ==="
base=$(run_gate)
printf '  %-58s %s\n' "unmutated tree" "$base"
[ "$base" = GREEN ] || { echo "baseline is not green; fix that before mutating"; exit 1; }

echo "=== mutations ==="
C=crates/contracts/operations
P=crates/providers

# 1. Contract moved to a domain its name no longer matches.
swap "contract name disagrees with its domain" \
    "$C/mesh-section/Cargo.toml" \
    'domain = "mesh.section"' 'domain = "mesh.slice"' RED

# 2. Provider moved to a domain its prefix no longer matches.
swap "provider name disagrees with its domain" \
    "$P/pointcloud/sdf/Cargo.toml" \
    'domain = "pointcloud.reconstruction"' 'domain = "pointcloud.meshing"' RED

# 3. A contract re-roled as a provider: it now carries `-contract` without
#    being one, and its bare name reservation disappears.
swap "contract demoted to a provider role" \
    "$C/mesh-boolean/Cargo.toml" \
    'role = "contract.operation"' 'role = "provider.mesh"' RED

# 4. A provider promoted to contract role while keeping an engine suffix.
swap "provider promoted to a contract role" \
    "$P/mesh/boolmesh/Cargo.toml" \
    'role = "provider.mesh"' 'role = "contract.operation"' RED

# 5. DECOY: a comment appended to the domain value. Looks like a metadata
#    edit, changes no meaning -- must stay green.
swap "comment appended to a domain (must NOT trip)" \
    "$P/mesh/boolmesh/Cargo.toml" \
    'domain = "mesh.boolean"' 'domain = "mesh.boolean" # decoy' GREEN

# 6. DECOY: the grandfathered orchestrator. If this trips, the exemption
#    list broke and the gate would demand a breaking rename of a crate
#    already published to crates.io.
swap "grandfathered orchestrator (must NOT trip)" \
    "crates/execution/compile/Cargo.toml" \
    'domain = "graph.compile"' 'domain = "graph.compile" # decoy' GREEN

echo "=== restored ==="
printf '  %-58s %s\n' "tree after restore" "$(run_gate)"

# Verify THIS SCRIPT reverted its own mutations, not that the tree is
# pristine: unrelated work in progress is the normal case when running this.
if [ "$(snapshot)" = "$BEFORE" ]; then
    echo "  manifests identical to pre-probe state"
else
    echo "  DIRTY -- restore failed; probe mutations leaked:"
    diff <(printf '%s\n' "$BEFORE") <(snapshot) | sed 's/^/    /'
    fail=1
fi

echo
[ "$fail" -eq 0 ] && echo "NAMING MATRIX PASSED" || echo "NAMING MATRIX FAILED"
exit "$fail"
