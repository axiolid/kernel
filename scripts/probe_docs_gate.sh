#!/usr/bin/env bash
# Mutation probe: can `cargo xtask docs --check` actually fail?
# Each mutation makes a generated page stale or breaks a README rule; the
# check must reject every one, for the right reason, and accept the
# untouched tree.
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
cargo build --quiet -p xtask || exit 1
fail=0
created=()
BAK="${TMPDIR:-/tmp}/docsprobe.$$"
mkdir -p "$BAK"
cleanup() {
  for f in "${created[@]}"; do rm -f "$f"; done
  for f in "$BAK"/*.path; do [ -e "$f" ] && cp "${f%.path}" "$(cat "$f")"; done
  rm -rf "$BAK"
}
trap cleanup EXIT
check() {
  printf "  %-58s" "$1"
  local out got
  out=$(cargo xtask docs --check 2>&1)
  got=$?
  if [ "$2" = "$got" ] && { [ -z "${3:-}" ] || grep -qF -- "$3" <<<"$out"; }; then
    echo ok
  else
    echo "MISS (want=$2${3:+ \"$3\"} got=$got)"; fail=1
  fi
  cleanup
  created=()
  mkdir -p "$BAK"
}
save() { local key; key=$(echo "$1" | tr '/' '_'); cp "$1" "$BAK/$key"; echo "$1" >"$BAK/$key.path"; }
new() { mkdir -p "$(dirname "$1")"; printf '%s\n' "$2" >"$1"; created+=("$1"); }

CRATE=crates/foundation/core
PAGE=docs/reference/crates/axiolid-core.md
REF='- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-core)'
echo "=== baseline ==="
check "untouched tree passes" 0
echo "=== generated pages ==="
save "$PAGE"; printf 'hand edit\n' >>"$PAGE"
check "hand-edited crate page" 1 "$PAGE is stale"
save docs/reference/index.md; rm docs/reference/index.md
check "missing index" 1 "docs/reference/index.md is stale"
new docs/reference/crates/axiolid-gone.md "# axiolid-gone"
check "page for a crate that no longer exists" 1 "axiolid-gone.md is no longer generated"
save docs/.vitepress/data/facts.json; printf '{}\n' >docs/.vitepress/data/facts.json
check "stale sidebar facts" 1 "facts.json is stale"
save docs/architecture/crate-map.md; printf 'x\n' >>docs/architecture/crate-map.md
check "stale architecture crate map" 1 "crate-map.md is stale"
save "$CRATE/README.md"; printf '\n## Design notes\n\nA new note.\n' >>"$CRATE/README.md"
check "README change not regenerated" 1 "$PAGE is stale"
save "$CRATE/CHANGELOG.md"; printf '\n## [9.9.9] - 2099-01-01\n\n- probe\n' >>"$CRATE/CHANGELOG.md"
check "release not regenerated (page)" 1 "$PAGE is stale"
save "$CRATE/CHANGELOG.md"; printf '\n## [9.9.9] - 2099-01-01\n\n- probe\n' >>"$CRATE/CHANGELOG.md"
check "release not regenerated (changelog page)" 1 "docs/reference/changelog.md is stale"
save "$CRATE/CHANGELOG.md"; sed -i 's/^## \[Unreleased\]$/&\n\n### Added\n\n- pending/' "$CRATE/CHANGELOG.md"
check "an [Unreleased] entry needs no regeneration" 0
echo "=== READMEs ==="
save "$CRATE/README.md"; grep -vF -- "$REF" "$CRATE/README.md" >"$BAK/r" && cp "$BAK/r" "$CRATE/README.md"
check "README without its reference link" 1 "missing the line \`$REF\`"
save "$CRATE/README.md"; sed -i 's#reference/crates/axiolid-core)#reference/crates/axiolid-mesh)#' "$CRATE/README.md"
check "README linking another crate's page" 1 "missing the line \`$REF\`"
save "$CRATE/README.md"; sed -i 's#docs.rs/axiolid-core)#docs.rs/axiolid)#' "$CRATE/README.md"
check "README without its docs.rs link" 1 "docs.rs/axiolid-core"
save "$CRATE/README.md"; printf 'This crate is not yet published.\n' >>"$CRATE/README.md"
check "README publication claim" 1 "publication claim"
save "$CRATE/README.md"; printf 'Stable since 0.3.1.\n' >>"$CRATE/README.md"
check "README version claim" 1 "version claim"
[ "$fail" = 0 ] && echo "DOCS_PROBE=PASS" || echo "DOCS_PROBE=FAIL"
exit "$fail"
