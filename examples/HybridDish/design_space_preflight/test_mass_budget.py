#!/usr/bin/env python3
"""Mass-budget unit tests A-E. Importable and runnable before revalidation."""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import run_transport_screen as rts  # noqa: E402


if __name__ == "__main__":
    rts.verify_random_drive()
    rts.run_mass_budget_unit_tests()
