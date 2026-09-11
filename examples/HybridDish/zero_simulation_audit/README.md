# Zero-simulation evidence audit

Analysis-only package for frozen NARMA, waveform, BenchA and BenchA2 data.
It invokes no Java or BSim and writes only to this directory plus the sibling
`external_analysis/` provenance archive.

Run from the repository root:

```powershell
python examples/HybridDish/zero_simulation_audit/run_zero_simulation_audit.py
```

The command reruns Stage 0 sanity, all eligible analyses, CSV generation,
the IPC affordability note, provenance copying/hash verification, and the
integrated report. See [PROTOCOL.md](./PROTOCOL.md) for the rules frozen before
calculation and [ZERO_SIMULATION_AUDIT.md](./ZERO_SIMULATION_AUDIT.md) for the
findings.
