#!/usr/bin/env bash
set -euo pipefail

cargo +1.88.0 xtask ffi check
cargo +1.88.0 build --release -p axiolid-capi

target_dir="${CARGO_TARGET_DIR:-target}"
output="${target_dir}/axiolid-capi-smoke"
if [[ "$(uname -s)" == Linux ]]; then
  header_symbols="$(grep -oE 'axiolid_v0_4_[A-Za-z0-9_]+' crates/facade/axiolid-capi/include/axiolid.h | sort -u)"
  library_symbols="$(nm -D --defined-only "${target_dir}/release/libaxiolid_capi.so" | awk '$2 == "T" {print $3}' | sort -u)"
  diff -u <(printf '%s\n' "$header_symbols") <(printf '%s\n' "$library_symbols")

  # The cdylib must carry a SONAME. Without one a consumer's DT_NEEDED entry
  # records whatever path the linker saw -- typically relative -- and ELF
  # never resolves a NEEDED value containing '/' through RPATH/RUNPATH, so
  # the consumer only loads from the build tree (kernel#70).
  soname="$(readelf -d "${target_dir}/release/libaxiolid_capi.so" \
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
