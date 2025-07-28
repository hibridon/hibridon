# AIRY

The flag `AIRYFL` should be set to 
- .TRUE. if the Airy integrator is to be used
- .FALSE. if the Airy integrator is not to be used

If `AIRYFL` is .FALSE., then the parameters which govern the Airy integration (`AIRYFL`, [`TOLAI`](./tolai.md), [`RINCR`](./rincr.md) and [`RENDAI`](./rendai.md) are assigned the following values:
- `PRAIRY` = .FALSE.
- `TOLAI` = 1.
- `RINCR` = 1.
- `RENDAI` = `RENDLD`

Note that the last assignment insures that if `AIRYFL` = .FALSE., then the [`CC`](../Close-coupled-equations.md) equations are integrated from [`RSTART`](./rstart.md) to [`RENDLD`](./rendld.md) only.

# PRAIRY

If the flag `PRAIRY` is set .TRUE., then at each step of the Airy integration the values of RNOW (the mid-point of the current interval), `DRNOW` (the width of the current interval), `CDIAG` (the maximum estimated diagonal correction), and `COFF` (the maximum estimated off-diagonal correction) are printed out. The definition of `CDIAG` and `COFF` is explained in M. H. Alexander, J. Chem. Phys. **81**, 4510 (1984).

⚠️ Do not set `AIRYFL` and [`LOGDFL`](./logdfl.md) both .FALSE.

See also the [`RINCR`](./rincr.md), [`LOGDFL`](./logdfl.md) and [`TOLAI`](./tolai.md) help files.