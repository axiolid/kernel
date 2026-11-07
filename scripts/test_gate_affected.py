"""Unit tests for scripts/gate-affected.py (#443): which paths are gate
infrastructure, and which workspace packages a change affects."""
import importlib.util
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "gate-affected.py"
spec = importlib.util.spec_from_file_location("gate_affected", SCRIPT)
gate_affected = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gate_affected)


def package(name, directory, deps=()):
    return {
        "name": name,
        "manifest_path": f"/ws/{directory}/Cargo.toml",
        "dependencies": [{"name": dep, "kind": None} for dep in deps],
    }


# base <- mid <- top <- app;  base <- side;  lone;  app also needs external `serde`.
# Dependents are listed first, so one pass over the list cannot find a
# transitive dependent by luck.
METADATA = {
    "workspace_root": "/ws",
    "packages": [
        package("app", "crates/facade/app", ["top", "serde"]),
        package("top", "crates/algorithms/top", ["mid"]),
        package("mid", "crates/contracts/mid", ["base"]),
        package("side", "crates/algorithms/side", ["base"]),
        package("base", "crates/foundation/base"),
        package("lone", "crates/algorithms/lone"),
        package("nested", "crates/algorithms/lone/nested", ["lone"]),
    ],
}


class InfraTests(unittest.TestCase):
    def test_infrastructure_paths(self):
        for path in (
            "scripts/gate.sh",
            "tools/xtask/src/main.rs",
            "Cargo.toml",
            "Cargo.lock",
            "crates/algorithms/top/Cargo.toml",
            "crates/facade/axiolid-capi/src/lib.rs",
            "native/CMakeLists.txt",
            "tests/consumers/rust-facade-application/src/main.rs",
            "tests/downstream/test_downstream_consumers.py",
            ".github/workflows/ci.yml",
            "docs/architecture/capability-ledger.toml",
        ):
            self.assertTrue(gate_affected.is_infra(path), path)

    def test_ordinary_paths_are_not_infrastructure(self):
        for path in (
            "crates/algorithms/top/src/lib.rs",
            "crates/algorithms/top/README.md",
            "docs/guide/contributing.md",
            "docs/architecture/crate-map.md",
            "tests/native/test_native_packaging.py",
            "tools/benchmark/src/main.rs",
        ):
            self.assertFalse(gate_affected.is_infra(path), path)

    def test_infra_files_filters_in_order(self):
        files = ["a/b.rs", "scripts/x.sh", "Cargo.lock"]
        self.assertEqual(gate_affected.infra_files(files), ["scripts/x.sh", "Cargo.lock"])


class AffectedTests(unittest.TestCase):
    def affected(self, *files):
        return gate_affected.affected_packages(METADATA, list(files))

    def test_leaf_change_affects_only_itself(self):
        self.assertEqual(self.affected("crates/facade/app/src/lib.rs"), ["app"])

    def test_every_transitive_reverse_dependent_is_included(self):
        self.assertEqual(
            self.affected("crates/foundation/base/src/lib.rs"),
            ["app", "base", "mid", "side", "top"],
        )

    def test_middle_change_pulls_in_dependents_not_dependencies(self):
        self.assertEqual(self.affected("crates/contracts/mid/src/a.rs"), ["app", "mid", "top"])

    def test_nested_package_owns_its_files(self):
        self.assertEqual(self.affected("crates/algorithms/lone/nested/src/lib.rs"), ["nested"])
        self.assertEqual(self.affected("crates/algorithms/lone/src/lib.rs"), ["lone", "nested"])

    def test_files_outside_every_package_affect_nothing(self):
        self.assertEqual(self.affected("docs/guide/x.md", "README.md"), [])

    def test_directory_prefix_is_not_a_substring_match(self):
        self.assertEqual(self.affected("crates/algorithms/lonely/src/lib.rs"), [])

    def test_external_dependencies_are_ignored(self):
        self.assertEqual(self.affected("serde/src/lib.rs"), [])


class CommandLineTests(unittest.TestCase):
    def run_main(self, command, files):
        with tempfile.NamedTemporaryFile("w", suffix=".txt") as listing:
            listing.write("\n".join(files))
            listing.flush()
            return gate_affected.main([command, "--files-from", listing.name])

    def test_infra_exit_codes(self):
        self.assertEqual(self.run_main("infra", ["scripts/gate.sh"]), 0)
        self.assertEqual(self.run_main("infra", ["crates/algorithms/top/src/lib.rs"]), 1)

    def test_unknown_base_is_an_error_not_an_empty_answer(self):
        self.assertEqual(gate_affected.main(["infra", "--base", "no-such-ref-443"]), 2)


if __name__ == "__main__":
    unittest.main()
