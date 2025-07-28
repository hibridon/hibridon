The command MRCRS allows you to calculate the m-state resolved cross sections defined by Eqs. (14) and (28) of B. Follmeg, P. Rosmus, and H.-J. Werner, J. Chem. Phys. **93**, 4687 (1990). The command line is
```
MRCRS,JOB,IENERG
```
where the input parameters are:
  - JOB,IENERG  The program searches for previous calculated tensor cross sections in the file Job{ienerg}.tcb, where IENERG is an integer (IENERG = 1, 2, ... etc.)

The MRCRS command must have been preceded by the command [`TENXSC`]../param/tenxsc.md), which creates the table of tensor cross sections used.

⚠️ Note that the m-state resolved cross sections for a transition from rotational level J to J' will be correct only if the parameter KMAX in the TENXSC command is
```
KMAX .ge. J+J'
```
MRCRS invokes the calculation of the m-state resolved cross sections appropriate to a Legendre expansion of a cylindrically symmetric distribution of the relative velocity vectors with respect to the z-axis. Here too, the calculated Legendre terms will be correct only for values of $\chi$ (the Legendre index) which are .le. the parameter LAMMAX in the previous invocation of TENXSC.
The calculated m-state resolved cross sections are printed out to the [file](../Files.md) Job{ienerg}.mcs as well as to the standard output