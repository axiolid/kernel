from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


pack = load("package_native", ROOT / "scripts/package-native.py")
verify = load("verify_native", ROOT / "scripts/verify-native-package.py")


class PackageBudgetTests(unittest.TestCase):
    def test_packager_budget_equals_verifier_budget(self) -> None:
        # The packager checks the verifier's limits early; a copy that drifts
        # would pass packaging and fail verification again (#431).
        self.assertEqual(pack.MAX_MEMBER_BYTES, verify.MAX_MEMBER_BYTES)
        self.assertEqual(pack.MAX_ARCHIVE_BYTES, verify.MAX_ARCHIVE_BYTES)


if __name__ == "__main__":
    unittest.main()
