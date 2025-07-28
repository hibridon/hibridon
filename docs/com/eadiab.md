# EADIAB

The command `EADIAB` requests the determination of the adiabatic energies as a function of $R$. These are defined by Eq. (6) of the help link which summarizes the [close coupling](../Close-coupled-equations.md) method.
The command line syntax can be one of the following:
```
EADIAB
EADIAB,<jobname>
EADIAB,<jobname>,<nch_max>
EADIAB,<jobname>,<nch_min>,<nch_max>
```
where
- `jobname`:  the jobname under which the wavefunction information from the previous run is stored, as the file `<jobname>.wfu`. If missing, the current `jobname` will be used. The default `jobname` in Hibridon is `job`
- `nch_min`: the ID of the first adiabatic energy for each $R$ to be printed, assuming the adiabatic energies are labeled starting from the lowest one. If missing, the default value 1 will be used.
- `nch_max`: the ID of the last adiabatic energy for each $R$ to be printed. If missing, the default value 10 will be used. Specifically, if `nch_max` is set to 0, all adiabatic energies will be printed.

The adiabatic energies will be stored in file `<jobname>.eadiab`, which is a text file with one header line. Units of $R$ and adiabatic energies used in the file are Bohr and wavenumber, respectively.

⚠️ In order to create the `<jobname>.wfu` file, the previous run must have been done with the flags [`WAVEFL`](../param/wavefl.md) and [`AIRYFL`](../param/airyfl) both set `.TRUE.`. In that run, [`WRSMAT`](../param/wrsmat.md) controls whether information other than the adiabatic energies, *e.g.* the transformation matrices, will be written to the `wfu` file. If [`WRSMAT`](../param/wrsmat.md) is set to FALSE, the `wfu` file will be much smaller but cannot be used for [`PSI`](./psi.md) or [`FLUX`](./flux.md) calculations.

The only adiabatic energies which are written out are those which correlate adiabatically at large $R$ with the internal states of the collision partners whose rotational angular momenta and additional indices are specified in the arrays [`JOUT`](./jout.md) and [`INDOUT`](./indout.md). Specifically, adiabatic energies are given for each channel for which
```
J(i) = JOUT(i)
```
and for which the vector element `INDOUT(i)` corresponds to the quantum number L and the extra index of the desired channel packed into a single array as follows:
```
INDOUT(i) = sign[IND(i)] x [100 L(i) + | IND(i) | ]
```
The asymptotic energy ordering of the internal states can be determined by running with the flag [`BASTST`]../param/bastst.md) set .TRUE.

⚠️ The adiabatic energies are determined at all values of $R$ lying between [`RENDLD`]RSTART,-RENDLD,-and-RENDAI) and [`RENDAI`]RSTART,-RENDLD,-and-RENDAI), with the grid size and spacing controlled by the same parameters which govern the AIRY integration, namely [`SPAC`]../param/spac.md), [`FSTFAC`]../param/fstfac.md), [`TOLAI`]TOLAI-and-RINCR), and [`RINCR`]TOLAI-and-RINCR).

⚠️ Highly localized avoided crossings of the adiabatic curves often occur. Use of a fine grid in $R$ is the only way to interpret properly the adiabatic energies in these cases.