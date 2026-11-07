"""Tests for the tier, meta-probe and slot logic of scripts/gate.sh (#443).

Every run is a dry run (GATE_DRY=1: each step is a sleep) in a private gate
root and lock directory, so nothing builds, nothing takes a real slot, and a
gate running this test does not queue behind itself.
"""
import os
import signal
import subprocess
import tempfile
import time
import unittest
from pathlib import Path

GATE = Path(__file__).resolve().parent / "gate.sh"


def _alive(pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    return True


class GateCase(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.dir = Path(self.temp.name)

    def env(self, changed=(), **extra):
        listing = self.dir / "changed.txt"
        listing.write_text("\n".join(changed))
        # A gate running these tests exports its own slot's settings; the
        # gate under test must start from none of them.
        inherited = ("CARGO_BUILD_JOBS", "RUST_TEST_THREADS", "CARGO_TARGET_DIR")
        env = {
            key: value
            for key, value in os.environ.items()
            if not key.startswith("GATE_") and key not in inherited
        }
        env.update(
            GATE_DRY="1",
            GATE_ROOT=str(self.dir / "root"),
            GATE_LOCK_DIR=str(self.dir / "locks"),
            GATE_CHANGED_FILES_FROM=str(listing),
        )
        env.update(extra)
        return env

    def start(self, *args, **kwargs):
        return subprocess.Popen(
            [str(GATE), *args],
            env=self.env(**kwargs),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            start_new_session=True,
        )

    def run_gate(self, *args, **kwargs):
        process = self.start(*args, **kwargs)
        out, _ = process.communicate(timeout=120)
        return process.returncode, out


class TierTests(GateCase):
    def test_branch_tier_on_a_crate_local_change(self):
        code, out = self.run_gate("--branch", changed=["crates/algorithms/exact/src/lib.rs"])
        self.assertEqual(code, 0, out)
        self.assertIn("gate tier: branch\n", out)
        self.assertIn("build (affected crates)", out)
        self.assertNotIn("build --workspace", out)
        self.assertNotIn("architecture mutation probe", out)
        self.assertIn("gate meta probes: skipped", out)
        self.assertIn("GATE PASSED", out)

    def test_branch_tier_with_no_crate_changed_skips_the_cargo_steps(self):
        code, out = self.run_gate("--branch", changed=["docs/guide/contributing.md"])
        self.assertEqual(code, 0, out)
        self.assertIn("skip (no workspace crate changed)", out)
        self.assertNotIn("build (affected crates)", out)

    def test_branch_tier_escalates_when_gate_infrastructure_changed(self):
        for path in (
            "scripts/gate.sh",
            "tools/xtask/src/main.rs",
            "crates/algorithms/exact/Cargo.toml",
            "native/CMakeLists.txt",
            "docs/architecture/closure-profiles.toml",
        ):
            with self.subTest(path=path):
                code, out = self.run_gate("--branch", changed=[path])
                self.assertEqual(code, 0, out)
                self.assertIn("escalated from --branch", out)
                self.assertIn("build --workspace", out)
                self.assertNotIn("build (affected crates)", out)

    def test_full_gate_skips_meta_probes_without_infrastructure_change(self):
        code, out = self.run_gate(changed=["crates/algorithms/exact/src/lib.rs"])
        self.assertEqual(code, 0, out)
        self.assertIn("gate tier: full\n", out)
        self.assertIn("build --workspace", out)
        self.assertIn("architecture mutation probe", out)
        self.assertIn("skip (no gate infrastructure change", out)
        self.assertIn("gate meta probes: skipped", out)

    def test_full_gate_runs_meta_probes_when_infrastructure_changed(self):
        code, out = self.run_gate(changed=["scripts/gate-lib.sh"])
        self.assertEqual(code, 0, out)
        self.assertIn("gate meta probes: run", out)
        self.assertNotIn("skip (", out)

    def test_gate_meta_forces_the_probes(self):
        code, out = self.run_gate(changed=[], GATE_META="1")
        self.assertEqual(code, 0, out)
        self.assertIn("gate meta probes: run", out)
        self.assertNotIn("skip (", out)

    def test_undecidable_change_counts_as_infrastructure_changed(self):
        code, out = self.run_gate(
            "--branch", GATE_CHANGED_FILES_FROM=str(self.dir / "does-not-exist")
        )
        self.assertEqual(code, 0, out)
        self.assertIn("escalated from --branch", out)
        self.assertIn("gate meta probes: run", out)

    def test_unknown_argument_is_refused(self):
        code, out = self.run_gate("--bogus")
        self.assertEqual(code, 2, out)

    def test_failure_keeps_the_exit_code_and_the_verdict(self):
        code, out = self.run_gate(GATE_DRY_FAIL="fmt --check")
        self.assertEqual(code, 1, out)
        self.assertIn("GATE FAILED", out)
        self.assertNotIn("GATE PASSED", out)

    def test_the_log_file_holds_the_output(self):
        log = self.dir / "own.log"
        code, out = self.run_gate(GATE_LOG=str(log))
        self.assertEqual(code, 0, out)
        self.assertIn(f"gate log: {log}", out)
        self.assertIn("GATE PASSED", log.read_text())


class SlotTests(GateCase):
    SLOW = dict(GATE_SLOTS="1", GATE_DRY_SLEEP="0.15", GATE_SLOT_POLL="0.1")

    def wait_for(self, process, text, seconds=30):
        deadline = time.time() + seconds
        seen = ""
        while time.time() < deadline:
            line = process.stdout.readline()
            if not line:
                break
            seen += line
            if text in line:
                return seen
        self.fail(f"never saw {text!r}:\n{seen}")

    def test_each_slot_builds_in_its_own_target_directory(self):
        code, out = self.run_gate()
        self.assertEqual(code, 0, out)
        self.assertIn(f"CARGO_TARGET_DIR={self.dir / 'root' / 'slot-1'})", out)

    def test_a_run_waits_while_every_slot_is_taken(self):
        first = self.start(**self.SLOW)
        self.addCleanup(first.wait)
        self.wait_for(first, "gate slot: 1 of 1")
        second = self.start(**self.SLOW)
        out, _ = second.communicate(timeout=120)
        first_out, _ = first.communicate(timeout=120)
        self.assertEqual(second.returncode, 0, out)
        self.assertIn("all 1 taken", out)
        self.assertIn("waiting", out)
        self.assertIn("gate slot: 1 of 1", out)
        self.assertNotIn("all 1 taken", first_out)
        self.assertIn("GATE PASSED", first_out)

    def test_two_slots_run_two_gates_and_queue_the_third(self):
        slow = dict(self.SLOW, GATE_SLOTS="2")
        first = self.start(**slow)
        self.addCleanup(first.wait)
        self.wait_for(first, "gate slot: 1 of 2")
        second = self.start(**slow)
        self.addCleanup(second.wait)
        self.wait_for(second, "gate slot: 2 of 2")
        third = self.start(**slow)
        out, _ = third.communicate(timeout=120)
        first.communicate(timeout=120)
        second.communicate(timeout=120)
        self.assertIn("all 2 taken", out)
        self.assertEqual(third.returncode, 0, out)

    def test_threads_are_split_between_slots_unless_preset(self):
        share = max(1, (os.cpu_count() or 1) // 2)
        _, out = self.run_gate(GATE_SLOTS="2")
        self.assertIn(f"CARGO_BUILD_JOBS={share} RUST_TEST_THREADS={share}", out)
        _, out = self.run_gate(GATE_SLOTS="2", CARGO_BUILD_JOBS="7", RUST_TEST_THREADS="5")
        self.assertIn("CARGO_BUILD_JOBS=7 RUST_TEST_THREADS=5", out)

    def test_steps_run_niced_unless_disabled(self):
        _, out = self.run_gate()
        self.assertIn("nice -n 10", out)
        _, out = self.run_gate(GATE_NICE="0")
        self.assertIn("unniced", out)

    def test_a_failed_gate_releases_its_slot(self):
        code, out = self.run_gate(GATE_SLOTS="1", GATE_DRY_FAIL="fmt --check")
        self.assertEqual(code, 1, out)
        code, out = self.run_gate(GATE_SLOTS="1")
        self.assertEqual(code, 0, out)
        self.assertNotIn("taken", out)

    def test_a_killed_gate_releases_its_slot(self):
        first = self.start(**self.SLOW)
        self.wait_for(first, "gate slot: 1 of 1")
        os.killpg(first.pid, signal.SIGTERM)
        first.communicate(timeout=60)
        code, out = self.run_gate(GATE_SLOTS="1")
        self.assertEqual(code, 0, out)
        self.assertNotIn("taken", out)

    def test_a_straggling_step_process_does_not_keep_the_slot(self):
        first = self.start(GATE_SLOTS="1", GATE_DRY_SLEEP="20")
        self.wait_for(first, "gate slot: 1 of 1")
        time.sleep(1)
        listing = subprocess.run(
            ["ps", "-o", "pid=,comm=", "-g", str(first.pid)],
            capture_output=True,
            text=True,
        ).stdout.split("\n")
        rows = [line.split(None, 1) for line in listing if line.strip()]
        stragglers = [int(pid) for pid, comm in rows if comm.strip() == "sleep"]
        self.assertTrue(stragglers, listing)
        self.addCleanup(lambda: [os.kill(pid, signal.SIGKILL) for pid in stragglers if _alive(pid)])
        # The gate dies; the step's own process (a dry-run sleep) lives on.
        for pid, comm in rows:
            if comm.strip() != "sleep":
                os.kill(int(pid), signal.SIGKILL)
        first.communicate(timeout=60)
        code, out = self.run_gate(GATE_SLOTS="1", GATE_DRY_SLEEP="0")
        self.assertEqual(code, 0, out)
        self.assertNotIn("taken", out)

    def test_a_slot_last_built_from_another_tree_rebuilds_this_one(self):
        slot = self.dir / "root" / "slot-1"
        slot.mkdir(parents=True)
        (slot / ".gate-tree").write_text("/elsewhere\n")
        code, out = self.run_gate(GATE_SLOTS="1")
        self.assertEqual(code, 0, out)
        self.assertIn("last built from /elsewhere, rebuilding this tree's", out)
        # The slot now belongs to this tree: the next gate here reuses it.
        code, out = self.run_gate(GATE_SLOTS="1")
        self.assertEqual(code, 0, out)
        self.assertNotIn("rebuilding", out)

    def test_a_new_slot_rebuilds_this_tree(self):
        code, out = self.run_gate(GATE_SLOTS="1")
        self.assertEqual(code, 0, out)
        self.assertIn("last built from an unknown tree, rebuilding", out)

    def test_bad_slot_count_is_refused(self):
        for bad in ("0", "x", "-1"):
            with self.subTest(slots=bad):
                code, out = self.run_gate(GATE_SLOTS=bad)
                self.assertEqual(code, 1, out)
                self.assertIn("GATE_SLOTS must be a positive integer", out)


if __name__ == "__main__":
    unittest.main()
