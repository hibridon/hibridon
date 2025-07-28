# PSI

The command PSI allows you to determine the scattering wavefunction from data stored during a previous calculation. Specifically this command will determine and write out one column of the F(R) matrix defined in the help file which describes the [close coupling](https://www2.chem.umd.edu/groups/alexander/hibridon/hib43/closecoupled.html) method.
The command line syntax is line is
```
PSI,{jobname},,INITIAL-J,INITIAL-L,INITIAL-IND
```
where
- {jobname}:  the jobname under which the wavefunction information from the previous calculation is stored, as the file {jobname}.wfu   In order to create the {jobname}.wfu file, the previous calculation must have been carried out with the flags [`WAVEFL`]../param/wavefl.md) and [`AIRYFL`]../param/airyfl.md) both set .true.

- INITIAL-J:  the rotational angular momentum of the incoming channel

- INITIAL-L:  the orbital angular momentum of the incoming channel

- INITIAL-IND:  the value of additional index for the incoming channel

The INITIAL-J, INITIAL-L, and INITIAL-IND indices specify which column of the F(R) matrix is selected. In addition, the only components (row elements) of this column which are written out are those for which the rotational angular momenta and additional indices are specified in the arrays [`JOUT`]../com/jout.md) and [`INDOUT`]../param/indout.md)

All selected components of the scattering wavefunctions are written to the file {jobname}.psi. The output is two successive matrices, corresponding, respectively, to the real and imaginary parts of $\mathbf{F}(R)$. The rows correspond to the successive values of R and the columns, to the components of $\mathbf{F}(R)$ which have been selected as described in the preceding paragraph.

 ⚠️ The wavefunction is determined at all values of R lying between [`RSTART`]RSTART,-RENDLD,-and-RENDAI) and [`RENDAI`]RSTART,-RENDLD,-and-RENDAI), with the grid size and spacing controlled by the same parameters which govern the LOGD and AIRY integration, namely [`SPAC`]../param/spac.md), [`FSTFAC`]../param/fstfac.md), [`TOLAI`]../param/tolai.md), and [`RINCR`]../param/rincr.md). In the output file {jobname}.psithe points that correspond to the AIRY integration range (RENDLD < R < RENDAI) are designated with negative value of R.

 ⚠️ It is essential that [`SPAC`]../param/spac.md) be signficantly less than the DeBroglie wavelength (which can be determined with the command [`DEBROGLI`]../com/debrogli.md)) is order to get a meaningful representation of the wavefunction.