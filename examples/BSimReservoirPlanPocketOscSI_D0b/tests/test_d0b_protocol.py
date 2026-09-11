#!/usr/bin/env python3
"""Assert D0b protocol is frozen before traces. No identity-grid run."""

from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))


def _load(name: str) -> dict:
    with (HERE / "configs" / name).open(encoding="utf-8") as handle:
        return json.load(handle)


class TestD0bProtocolFrozen(unittest.TestCase):
    def test_protocol_markdown_exists_and_frozen(self) -> None:
        path = HERE / "PROTOCOL.md"
        self.assertTrue(path.exists())
        text = path.read_text(encoding="utf-8")
        self.assertIn("frozen_before_traces:** true", text)
        self.assertIn("AHL_KICK_005", text)
        self.assertIn("PREDECLARED_AHL_KICK", text)
        self.assertIn("FAIL_NO_IDENTITY", text)

    def test_json_frozen_flags(self) -> None:
        for name in ("protocol.json", "ensemble.json", "mu_grid.json"):
            data = _load(name)
            self.assertIs(data["frozen_before_traces"], True, msg=name)

    def test_primary_is_ahl_kick_005(self) -> None:
        ens = _load("ensemble.json")
        self.assertEqual(ens["label"], "PREDECLARED_AHL_KICK")
        self.assertEqual(ens["primary_id"], "AHL_KICK_005")
        primary = next(a for a in ens["arms"] if a["id"] == "AHL_KICK_005")
        self.assertEqual(primary["A"], 0.0)
        self.assertEqual(primary["I"], 0.0)
        self.assertEqual(primary["H_i"], 0.05)
        self.assertEqual(primary["H_e"], 0.05)
        self.assertEqual(primary["H_i_history"], 0.05)

    def test_required_controls_present(self) -> None:
        ens = _load("ensemble.json")
        ids = {a["id"] for a in ens["arms"]}
        self.assertEqual(
            ids,
            {
                "AHL_KICK_005",
                "SI_BASAL_PERTURB",
                "AHL_KICK_001",
                "FAILED_PRIOR_IVP",
                "AHL_KICK_020",
            },
        )
        basal = next(a for a in ens["arms"] if a["id"] == "SI_BASAL_PERTURB")
        self.assertEqual(basal["I"], 1.0)
        self.assertEqual(basal["H_i"], 0.0)
        self.assertEqual(basal["H_i_history"], 0.0)

    def test_identity_mu_line_frozen(self) -> None:
        grid = _load("mu_grid.json")
        ident = [a["mu"] for a in grid["arms"] if a["role"] == "IDENTITY"]
        self.assertEqual(ident, [0.32, 0.36, 0.40, 0.44, 0.48, 0.52, 0.56, 0.60])
        extras = [a["mu"] for a in grid["arms"] if a["role"] == "FINITE_INTERVAL_EXTRA"]
        self.assertEqual(extras, [0.0, 0.10, 0.20, 0.28, 0.80, 1.00, 1.20, 1.50, 2.00])

    def test_peak_rules_copied_from_d0(self) -> None:
        p = _load("protocol.json")
        self.assertEqual(p["t_end"], 1000.0)
        self.assertEqual(p["t_discard"], 180.0)
        self.assertEqual(p["sample_dt"], 0.5)
        self.assertEqual(p["peak_min_spacing"], 25.0)
        self.assertEqual(p["peak_rel_prominence"], 0.20)
        self.assertEqual(p["min_peaks"], 4)
        self.assertEqual(p["identity_min_osc"], 6)
        self.assertEqual(p["spearman_min"], 0.95)

    def test_taken_constants_not_retuned(self) -> None:
        p = _load("params.json")
        self.assertEqual(p["alpha"], 2500.0)
        self.assertEqual(p["tau"], 10.0)
        self.assertEqual(p["gamma_H"], 0.01)
        self.assertEqual(p["k1"], 0.1)
        self.assertEqual(p["d"], 0.5)

    def test_d0_standing_not_rewritten_here(self) -> None:
        standing = HERE.parent / "PocketDish" / "D0_BULK_ORACLE_STANDING.md"
        self.assertTrue(standing.exists())
        text = standing.read_text(encoding="utf-8")
        self.assertIn("**Status: FAIL**", text)


if __name__ == "__main__":
    unittest.main()
