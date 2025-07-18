The command `EADIAB` requests the determination of the adiabatic energies as a function of R. These are defined by Eq. (6) of the help link which summarizes the [close coupling](Close-coupled-equations) method.
The command line syntax is
```
EADIAB,{jobname}
```
where
- {jobname}:  the jobname under which the wavefunction information from the previous run is stored, as the file {jobname}.wfu

⚠️ In order to create the {jobname}.wfu file, the previous run must have been done with the flags [WAVEFL](WAVEFL) and [AIRYFL](AIRYFL-and-PRAIRY) both set .TRUE.

The only adiabatic energies which are written out are those which correlate adiabatically at large R with the internal states of the collision partners whose rotational angular momenta and additional indices are specified in the arrays [JOUT](JOUT) and [INDOUT](INDOUT). Specifically, adiabatic energies are given for each channel for which
```
J(i) = JOUT(i)
```
and for which the vector element INDOUT(i) corresponds to the quantum number L and the extra index of the desired channel packed into a single array as follows:
```
INDOUT(i) = sign[IND(i)] x [100 L(i) + | IND(i) | ]
```
The asymptotic energy ordering of the internal states can be determined by running with the flag [BASTST](BASTST) set .TRUE.

⚠️ The adiabatic energies are determined at all values of R lying between [RENDLD](RSTART,-RENDLD,-and-RENDAI) and [RENDAI](RSTART,-RENDLD,-and-RENDAI), with the grid size and spacing controlled by the same parameters which govern the AIRY integration, namely [SPAC](SPAC), [FSTFAC](FSTFAC), [TOLAI](TOLAI-and-RINCR), and [RINCR](TOLAI-and-RINCR).

⚠️ Highly localized avoided crossings of the adiabatic curves often occur. Use of a fine grid in R is the only way to interpret properly the adiabatic energies in these cases.