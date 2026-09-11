# ChassisPocket

New named living dish that PocketDish / HybridDish can later import.
Not paper-2 HybridDish. Not a PocketDish Java graft.

I0 (locked PASS): locked Hertzian *E. coli* rods in a **closed**
PocketDish-A garage (`100 × 100 × 1 µm`). Standing:
[`CHASSISPOCKET_I0_STANDING.md`](CHASSISPOCKET_I0_STANDING.md).

I0b (scored **FAIL** — door not reached): same stack, T4/T5 particle
door open at frozen `W = 20 µm`. Standing:
[`CHASSISPOCKET_I0b_STANDING.md`](CHASSISPOCKET_I0b_STANDING.md).
Do not narrow `W` after seeing `N`. Protocol: [`PROTOCOL.md`](PROTOCOL.md).

```
compile_and_run.cmd smoke       ← I0 smoke N=32 (already scored)
compile_and_run.cmd             ← I0 full (already scored)
python check_i0.py

compile_and_run.cmd i0b smoke   ← I0b stop at N=32; filament check
compile_and_run.cmd i0b         ← I0b full 23400 s, seed 101
python check_i0b.py
```

AHL, Hill, NARMA, motility, death, and nutrient PDE are **OFF** in I0
and I0b.
