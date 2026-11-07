# shellcheck shell=bash
# Sourced by scripts/gate.sh (#443): tier choice, build slots, the step runner.
#
# Kept apart from gate.sh on purpose. bash reads a running script
# incrementally, so the mutation probes of this logic (probe_gate_tiers.py)
# must not edit gate.sh while a gate runs it; this file is read once, by
# `source`, before any step runs.

# Where slot build directories, locks' companions and logs live: next to the
# repository's common git dir, so every worktree of one clone shares them.
gate_root() {
  if [ -n "${GATE_ROOT:-}" ]; then printf '%s\n' "$GATE_ROOT"; return; fi
  local common
  common="$(git rev-parse --git-common-dir)" || return 1
  printf '%s/gate-target\n' "$(cd "$common/.." && pwd -P)"
}

# True when the change since the merge base touches gate infrastructure
# (scripts/gate-affected.py lists what counts). Exit 2 from it means the
# change cannot be determined; that counts as "changed", so a guess only ever
# runs more checks.
gate_infra_changed() {
  local rc=0
  infra_files="$(python3 scripts/gate-affected.py infra ${GATE_CHANGED_FILES_FROM:+--files-from "$GATE_CHANGED_FILES_FROM"} 2>&1)" || rc=$?
  case "$rc" in
    0) return 0 ;;
    1) return 1 ;;
    *) return 0 ;;
  esac
}

# Sets infra, infra_files, tier (branch|full) and meta (0|1). A --branch run
# escalates to the full gate when infrastructure changed; the gate's own
# mutation probes run on an infrastructure change or GATE_META=1.
gate_decide_tier() {
  infra=0
  infra_files=""
  gate_infra_changed && infra=1
  tier=full
  if [ "$1" -eq 1 ]; then
    tier=branch
    [ "$infra" -eq 1 ] && tier=full
  fi
  meta="$infra"
  [ "${GATE_META:-0}" = 1 ] && meta=1
}

# Take one of GATE_SLOTS (default 3) machine-wide slots with flock on fd 9 and
# build in that slot's own persistent target directory. The lock lives as long
# as this process: it is released on exit, failure and signal alike. Steps run
# with fd 9 closed (see `step`), so a straggling child cannot keep a slot.
gate_acquire_slot() {
  local slots="${GATE_SLOTS:-3}" n announced=0
  case "$slots" in '' | *[!0-9]* | 0)
    echo "GATE_SLOTS must be a positive integer, got '$slots'"
    return 1
    ;;
  esac
  mkdir -p "$GATE_LOCK_DIR" || return 1
  while :; do
    for n in $(seq 1 "$slots"); do
      exec 9>"$GATE_LOCK_DIR/slot-$n.lock" || return 1
      if flock -n 9; then
        GATE_SLOT="$n"
        CARGO_TARGET_DIR="$GATE_ROOT_DIR/slot-$n"
        export CARGO_TARGET_DIR
        mkdir -p "$CARGO_TARGET_DIR" || return 1
        echo "gate slot: $n of $slots (CARGO_TARGET_DIR=$CARGO_TARGET_DIR)"
        gate_claim_slot_tree "$CARGO_TARGET_DIR" || return 1
        gate_limit_resources
        return 0
      fi
    done
    if [ "$announced" -eq 0 ]; then
      echo "gate slot: all $slots taken (GATE_SLOTS), waiting for one to free up"
      announced=1
    fi
    sleep "${GATE_SLOT_POLL:-2}"
  done
}

# Make a slot build this worktree's sources, not another's. Cargo records a
# workspace member's sources by path relative to the workspace and judges
# them by modification time, so a slot last built from another worktree
# reuses that build for every file not modified here since: different code,
# and `env!("CARGO_MANIFEST_DIR")` naming the other tree. When the slot's
# last tree is not this one, every file of this worktree is touched, so each
# workspace member rebuilds from it; registry dependencies stay cached.
# GATE_DRY=1 only reports what it would touch.
gate_claim_slot_tree() {
  local owner="$1/.gate-tree" tree last=""
  tree="$(git rev-parse --show-toplevel)" || return 1
  [ -f "$owner" ] && last="$(cat "$owner")"
  if [ "$last" != "$tree" ]; then
    local count
    count="$(git -C "$tree" ls-files --cached --others --exclude-standard | wc -l)"
    echo "gate slot: last built from ${last:-an unknown tree}, rebuilding this tree's $count files"
    if [ -z "${GATE_DRY:-}" ]; then
      git -C "$tree" ls-files -z --cached --others --exclude-standard |
        (cd "$tree" && xargs -0 touch -c --) || return 1
    fi
  fi
  printf '%s\n' "$tree" >"$owner"
}

# Share the machine between slots: each gate gets nproc / GATE_SLOTS build and
# test threads (unless the caller set them) and runs below interactive work
# (nice 10, ionice idle where available). GATE_NICE=0 turns the niceness off.
gate_limit_resources() {
  local slots="${GATE_SLOTS:-3}" share
  share=$(($(nproc) / slots))
  [ "$share" -ge 1 ] || share=1
  export CARGO_BUILD_JOBS="${CARGO_BUILD_JOBS:-$share}"
  export RUST_TEST_THREADS="${RUST_TEST_THREADS:-$share}"
  gate_wrap=()
  if [ "${GATE_NICE:-1}" != 0 ]; then
    gate_wrap=(nice -n 10)
    command -v ionice >/dev/null 2>&1 && gate_wrap=(ionice -c3 nice -n 10)
  fi
  echo "gate resources: CARGO_BUILD_JOBS=$CARGO_BUILD_JOBS RUST_TEST_THREADS=$RUST_TEST_THREADS ${gate_wrap[*]:-unniced}"
}

# One gate step: print its name, run it, show the tail of its output on
# failure. GATE_DRY=1 replaces every command by a sleep (GATE_DRY_SLEEP
# seconds) so the queueing and tier logic can be tested without a build;
# GATE_DRY_FAIL=<name> makes the step of that exact name fail.
step() {
  local name="$1"
  shift
  printf '%-46s' "$name"
  if [ -n "${GATE_DRY:-}" ]; then
    set -- sleep "${GATE_DRY_SLEEP:-0}"
    [ "${GATE_DRY_FAIL:-}" = "$name" ] && set -- false
  fi
  local started=$SECONDS rc=0
  "${gate_wrap[@]}" "$@" 9>&- >"$gate_out" 2>&1 || rc=$?
  # Each step's wall time, so a slow step shows where the gate's time goes.
  if [ "$rc" -eq 0 ]; then
    echo "ok ($((SECONDS - started))s)"
  else
    echo "FAIL (exit $rc, $((SECONDS - started))s)"
    tail -25 "$gate_out" | sed 's/^/    /'
    fail=1
  fi
}

# A step that proves a gate itself can fail. It runs when gate infrastructure
# changed or GATE_META=1, otherwise it is skipped (the branch tier does not
# even mention it).
meta_step() {
  if [ "$meta" -eq 1 ]; then
    step "$@"
  elif [ "$tier" = full ]; then
    printf '%-46sskip (no gate infrastructure change; GATE_META=1 runs it)\n' "$1"
  fi
}
