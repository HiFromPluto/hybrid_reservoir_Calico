#!/usr/bin/env python3
"""Fast D0 equation and interpolant tests. No long traces. No NARMA."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parents[1]
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from danino_si_dde import (
    PRODUCTION_POCKETOSCSI_SPLIT,
    PRODUCTION_SI_HILL,
    DelayTape,
    DaninoSIDDE,
    load_json,
    production,
)


class TestD0Equations(unittest.TestCase):
    def test_si_hill_has_delta_in_numerator(self) -> None:
        # At large H the two forms agree; at H=0 both equal delta.
        self.assertAlmostEqual(
            production(0.0, 0.001, 2500.0, 0.1, PRODUCTION_SI_HILL), 0.001
        )
        # Distinct algebraic form: scale delta and compare.
        p_si = production(1.0, 10.0, 1.0, 1.0, PRODUCTION_SI_HILL)
        p_split = production(1.0, 10.0, 1.0, 1.0, PRODUCTION_POCKETOSCSI_SPLIT)
        self.assertAlmostEqual(p_si, (10.0 + 1.0) / 2.0)
        self.assertAlmostEqual(p_split, 10.0 + 1.0 / 2.0)
        self.assertNotAlmostEqual(p_si, p_split)

    def test_taken_p_values_are_close_but_form_differs(self) -> None:
        delta, alpha, k1 = 0.001, 2500.0, 0.1
        h = 1.0
        si = production(h, delta, alpha, k1, PRODUCTION_SI_HILL)
        split = production(h, delta, alpha, k1, PRODUCTION_POCKETOSCSI_SPLIT)
        self.assertGreater(abs(si - split), 0.0)
        self.assertLess(abs(si - split) / si, 1e-5)

    def test_delay_tape_raises_on_gap(self) -> None:
        tape = DelayTape(np.array([0.0, 0.0, 0.0, 0.0]), hi_history=1.0)
        with self.assertRaises(RuntimeError):
            tape.hi(1.0)

    def test_history_used_for_negative_time(self) -> None:
        tape = DelayTape(np.array([0.0, 1.0, 0.0, 0.0]), hi_history=0.25)
        self.assertEqual(tape.hi(-5.0), 0.25)
        self.assertEqual(tape.hi(0.0), 0.0)

    def test_params_match_si_table(self) -> None:
        p = load_json("params.json")
        self.assertEqual(p["C_A"], 1.0)
        self.assertEqual(p["C_I"], 4.0)
        self.assertEqual(p["delta"], 0.001)
        self.assertEqual(p["alpha"], 2500.0)
        self.assertEqual(p["tau"], 10.0)
        self.assertEqual(p["gamma_I"], 24.0)
        self.assertEqual(p["d"], 0.5)
        self.assertEqual(p["production_form"], "si_hill")

    def test_rhs_bulk_drops_d1(self) -> None:
        model = DaninoSIDDE()
        y = np.array([1.0, 2.0, 0.1, 0.2])
        dydt = model.rhs(y, h_tau=0.1, d=0.5, mu=1.5)
        self.assertEqual(dydt.shape, (4,))
        # He equation has no spatial term; changing Hi/He membrane part is finite.
        self.assertTrue(np.all(np.isfinite(dydt)))

    def test_narma_refused(self) -> None:
        from period_check import refuse_narma

        with self.assertRaises(SystemExit):
            refuse_narma(["run_d0.py", "--narma"])


if __name__ == "__main__":
    unittest.main()
