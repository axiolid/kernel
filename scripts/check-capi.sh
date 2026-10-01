#!/usr/bin/env bash
set -euo pipefail

cargo +1.88.0 xtask ffi check
cargo +1.88.0 build --release -p axiolid-capi

target_dir="${CARGO_TARGET_DIR:-target}"
output="${target_dir}/axiolid-capi-smoke"
header_symbols="$(grep -oE 'axiolid_v0_4_[A-Za-z0-9_]+' crates/facade/axiolid-capi/include/axiolid.h | sort -u)"

# kernel#56: this symbol-listing check used to run only under
# `uname -s == Linux`. The SONAME half of that block genuinely is an
# ELF-only concept, but the symbol-completeness diff got scoped under the
# same guard by accident, so on macOS and Windows the whole check -- the
# part that actually proves the cdylib exports every symbol the header
# declares -- silently never ran, and the script still exited 0. A gate
# with nothing to check must fail, not pass, so every branch below either
# performs a real diff or exits 1 with a reason; there is no silent no-op
# branch left.
os="$(uname -s)"
case "$os" in
  Linux)
    cdylib="${target_dir}/release/libaxiolid_capi.so"
    # GNU nm's dynamic symbol table: unmangled C names, no underscore
    # decoration. Restrict to defined text (function) symbols, same as
    # before this change.
    library_symbols="$(nm -D --defined-only "$cdylib" | awk '$2 == "T" {print $3}' | sort -u)"
    ;;
  Darwin)
    cdylib="${target_dir}/release/libaxiolid_capi.dylib"
    # BSD/Mach-O nm: `-g` restricts to external (global) symbols, `-U`
    # drops undefined ones, leaving only what the dylib actually defines
    # and exports. Mach-O additionally prefixes every C symbol with `_`
    # at the object-file level, which `nm` does not strip, so the leading
    # underscore has to be removed before comparing against the header.
    library_symbols="$(nm -gU "$cdylib" | grep -oE '_?axiolid_v0_4_[A-Za-z0-9_]+' | sed 's/^_//' | sort -u)"
    ;;
  *)
    # Covers Windows (MSYS/Git-Bash report MINGW64_NT-*/MSYS_NT-*) and any
    # other `uname -s` this script has not been taught yet. Try the
    # export-listing tools a Windows toolchain plausibly has on PATH
    # before giving up -- but giving up must still be a hard failure, not
    # a quiet pass, per kernel#56.
    cdylib="${target_dir}/release/axiolid_capi.dll"
    if command -v dumpbin >/dev/null 2>&1; then
      library_symbols="$(dumpbin /EXPORTS "$cdylib" | grep -oE 'axiolid_v0_4_[A-Za-z0-9_]+' | sort -u)"
    elif command -v llvm-nm >/dev/null 2>&1; then
      library_symbols="$(llvm-nm --defined-only "$cdylib" | grep -oE 'axiolid_v0_4_[A-Za-z0-9_]+' | sort -u)"
    else
      echo "check-capi.sh: no symbol-listing tool (nm/dumpbin/llvm-nm) found for 'uname -s' == '${os}'; cannot verify the cdylib exports every declared C ABI symbol. Install one of those tools or teach this script the platform instead of letting the check pass with nothing examined (kernel#56)." >&2
      exit 1
    fi
    ;;
esac

if [[ -z "$library_symbols" ]]; then
  echo "check-capi.sh: the symbol lister produced no defined axiolid_v0_4_* symbols from ${cdylib}; treating an empty result as a failed probe, not a pass (kernel#56)." >&2
  exit 1
fi
diff -u <(printf '%s\n' "$header_symbols") <(printf '%s\n' "$library_symbols")

# The cdylib must carry a SONAME. Without one a consumer's DT_NEEDED entry
# records whatever path the linker saw -- typically relative -- and ELF
# never resolves a NEEDED value containing '/' through RPATH/RUNPATH, so
# the consumer only loads from the build tree (kernel#70). SONAME is an
# ELF concept with no Mach-O/PE equivalent, so this stays Linux-only; it
# is a separate, narrower check from the symbol diff above and must not
# gate it again.
if [[ "$os" == Linux ]]; then
  soname="$(readelf -d "$cdylib" \
    | sed -n 's/.*SONAME.*\[\(.*\)\]/\1/p')"
  if [[ "$soname" != "libaxiolid_capi.so" ]]; then
    echo "capi cdylib SONAME is '${soname:-<none>}', expected libaxiolid_capi.so" >&2
    exit 1
  fi
fi
printf '#include "axiolid.h"\nint main() { return (int)AxiolidStatus_Ok; }\n' |
  "${CXX:-c++}" -std=c++17 -Wall -Wextra -Werror \
    -I crates/facade/axiolid-capi/include -x c++ -fsyntax-only -
"${CC:-cc}" -std=c11 -Wall -Wextra -Werror \
  -I crates/facade/axiolid-capi/include \
  crates/facade/axiolid-capi/tests/c/smoke.c \
  "${target_dir}/release/libaxiolid_capi.a" \
  -lm -ldl -lpthread -o "${output}"
"${output}"
