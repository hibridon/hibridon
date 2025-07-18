All files created by the Hibridon package are given lower case names with the first character in upper case.

All files written during a calculation, except the main output file, are designated {jobname}n.ext, where n denotes the energy at which the calculation has been carried out (n = 1 for first energy, n = 2 for second energy, etc.) and the file extension ext can be

  - `ics` for integral cross sections (created if the flag [WRXSEC](PRXSEC-and-WRXSEC) is set .true. and a necessary prerequisite to the command [PRINTC](PRINTC))

  - `xsc` for selected integral cross sections (created by the command [PRINTC](PRINTC))

  - `smt` for S matrix elements (created if the flag [WRSMAT](PRSMAT-and-WRSMAT) is set .true. or by the command [OPTIMIZE](OPTIMIZE); a necessary prerequisite to the commands [PRINTS](PRINTS) and [DIFCRS](DIFCRS))

  - `dcs` for differential cross sections (created by the command [DIFCRS](DIFCRS))
pcs for partial cross sections (created if the flag [WRPART](PRPART-and-WRPART) is set .true.)

  - `psc` for partial cross sections as output of the command [PARTC](PARTC))

  - `tcb` for a table of tensor cross sections (as output of the command [TENXSC](TENXSC) and necessary prerequisite for the command [MRCRS](MRCRS))

  - `tcs` for tensor cross sections (as output of the command [TENXSC](TENXSC))

  - `mcs` for a table of m-resolved integral cross sections (as output of the command [MRCRS](MRCRS))

  - `wfu` for wavefunction transformation matrices (created if the flag [WAVEFL](WAVEFL) is set .true. and a necessary prerequisite to the commands [FLUX](FLUX), [EADIAB](EADIAB), and [PSI](PSI))

  - `psi` for scattering wavefunctions (created by the command [PSI](PSI))

  - `flx` for fluxes or adiabatic energies (as output of the commands [FLUX](FLUX) and [EADIAB](EADIAB))

The jobname is set by the variable
```
JOB={jobname},
```
where {jobname} is a character string (not in quotes). The default value of {jobname} is Job

❗ Both {jobname} and the file extensions are converted internally to lower case with the first character of {jobname} in upper case.

❗ All of the above files are formatted, except the .smt, .tcb, and .wfu files which are unformatted (binary).