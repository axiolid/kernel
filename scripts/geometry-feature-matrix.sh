#!/usr/bin/env bash
# Feature-isolated geometry build matrix. Keep in sync with axiolid/Cargo.toml.
set -euo pipefail

cd "$(dirname "$0")/.."

step() {
    printf '%-52s' "$1"
    shift
    "$@" >/dev/null
    printf 'ok\n'
}

# A featureless facade must REFUSE to build (kernel#9). Asserting the failure,
# not tolerating it: if this ever succeeds, `default = []` has been undone or
# the diagnostic has been dropped, and consumers are silently paying for
# geometry again.
printf '%-46s' "axiolid facade: no features refuses"
if cargo check -q -p axiolid --no-default-features >/dev/null 2>&1; then
    printf 'FAIL (a featureless facade built; the compile_error is gone)\n'
    fail=1
else
    printf 'ok\n'
fi

# `standard` is the migration target the diagnostic names, so it has to work.
step "axiolid facade: standard bundle" cargo check -q -p axiolid --no-default-features --features standard

features=(
    mesh profiles curves surfaces topology primitives model
    evaluate nurbs tessellation spatial ray-mesh measure overlay field field-ops field-navigation heal
    contracts mesh-contracts mesh-boolean mesh-section graph-compile
    dispatch-mesh-boolean dispatch-mesh-section generate
    cpu parallel simd gpu
    standard discrete parametric advanced full
)
for feature in "${features[@]}"; do
    step "axiolid facade feature: ${feature}" \
        cargo check -q -p axiolid --no-default-features --features "$feature"
done

# `standard` is what a consumer gets by following the migration note, so its
# tests are the ones that stand in for "the default build" now that the real
# default compiles nothing.
step "axiolid facade: standard tests" cargo test -q -p axiolid --no-default-features --features standard
step "axiolid facade: all features" cargo test -q -p axiolid --all-features

step "common contracts" cargo test -q -p axiolid-contracts
step "guarantee vocabulary" cargo test -q -p axiolid-guarantees
step "shared mesh contracts" cargo test -q -p axiolid-mesh-contracts
step "tessellation contract" cargo test -q -p axiolid-tessellation-contract
step "mesh boolean contract" cargo test -q -p axiolid-mesh-boolean-contract
step "mesh section contract" cargo test -q -p axiolid-mesh-section-contract
step "mesh compile contract" cargo test -q -p axiolid-mesh-compile-contract
step "dispatch: identity only" cargo check -q -p axiolid-dispatch --no-default-features
step "dispatch: mesh boolean" cargo test -q -p axiolid-dispatch --no-default-features --features mesh-boolean
step "dispatch: mesh section" cargo test -q -p axiolid-dispatch --no-default-features --features mesh-section
step "dispatch: all" cargo test -q -p axiolid-dispatch --all-features

step "CPU context: portable" \
    cargo check -q -p axiolid-backend-cpu --no-default-features
step "CPU context: SIMD" \
    cargo check -q -p axiolid-backend-cpu --no-default-features --features simd
step "CPU context: parallel" \
    cargo check -q -p axiolid-backend-cpu --no-default-features --features parallel
step "CPU context: SIMD + parallel" \
    cargo test -q -p axiolid-backend-cpu --all-features
step "GPU adapter contract" cargo test -q -p axiolid-backend-gpu

if rustup target list --installed | grep -qx 'aarch64-linux-android'; then
    step "CPU context: AArch64 compile" \
        cargo check -q -p axiolid-backend-cpu --target aarch64-linux-android \
        --no-default-features --features simd
else
    printf '%-52s%s\n' "CPU context: AArch64 compile" "skip (target not installed)"
fi
