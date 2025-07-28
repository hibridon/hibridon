All files created by the Hibridon package are given lower case names with the first character in upper case.

All files written during a calculation, except the main output file, are designated {jobname}n.ext, where n denotes the energy at which the calculation has been carried out (n = 1 for first energy, n = 2 for second energy, etc.) and the file extension ext can be

  - `ics` for integral cross sections (created if the flag [`WRXSEC`](./param/wrxsec.md) is set .true. and a necessary prerequisite to the command [`PRINTC`](./com/printc.md))

  - `xsc` for selected integral cross sections (created by the command [`PRINTC`](./com/printc.md))

  - `smt` for S matrix elements (created if the flag [`WRSMAT`](./param/wrsmat.md) is set .true. or by the command [`OPTIMIZE`](./com/optimize.md); a necessary prerequisite to the commands [`PRINTS`](./com/prints.md) and [`DIFCRS`](./com/difcrs.md))

  - `dcs` for differential cross sections (created by the command [`DIFCRS`](./com/difcrs.md))
pcs for partial cross sections (created if the flag [`WRPART`](./param/wrpart.md) is set .true.)

  - `psc` for partial cross sections as output of the command [`PARTC`](./com/partc.md))

  - `tcb` for a table of tensor cross sections (as output of the command [`TENXSC`](./com/tenxsc.md) and necessary prerequisite for the command [`MRCRS`](./com/mrcrs.md))

  - `tcs` for tensor cross sections (as output of the command [`TENXSC`](./com/tenxsc.md))

  - `mcs` for a table of m-resolved integral cross sections (as output of the command [`MRCRS`](./com/mrcrs.md))

  - `wfu` for wavefunction transformation matrices (created if the flag [`WAVEFL`](./param/wavefl.md) is set .true. and a necessary prerequisite to the commands [`FLUX`](./com/flux.md), [`EADIAB`](./com/eadiab.md), and [`PSI`](./com/psi.md))

  - `psi` for scattering wavefunctions (created by the command [`PSI`](./com/psi.md))

  - `flx` for fluxes or adiabatic energies (as output of the commands [`FLUX`](./com/flux.md) and [`EADIAB`](./com/eadiab.md))

The jobname is set by the variable
```
JOB={jobname},
```
where {jobname} is a character string (not in quotes). The default value of {jobname} is Job

❗ Both {jobname} and the file extensions are converted internally to lower case with the first character of {jobname} in upper case.

❗ All of the above files are formatted, except the .smt, .tcb, and .wfu files which are unformatted (binary).