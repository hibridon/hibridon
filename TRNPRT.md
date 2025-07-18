The command TRNPRT allows you to calculate the effective cross sections $Q^n(E)$, required for the calculation of transport properties. See G. C. Maitland, M. Mustafa, W. A. Wakeham, and F. R. W. McCourt, Mol. Phys. **61**, 359 (1987); L. Monchick and S. Green, J. Chem. Phys. **63**, 2000 (1975); G. A. Parker and R. T Pack, J. Chem. Phys. **68**, 1585 (1978).
The command line is
```
TRNPRT,JOB,IENERG,IN1,IN2,JTOTMX,JMIN,JMAX
```
where the input parameters are:
- JOB,IENERG  The program searches for an S-matrix in the file Job{ienerg}.smt, where IENERG is an integer (IENERG = 1, 2, ... etc.)
The default values of JOB is JOB and of INERG is 1
- IN1,IN2  The initial and final values of the additional index
- JTOTMX  The maximum value included in the summation over total angular momentum (upper case $J$ in Eqs. 30 and 31 of the above article)
- JMIN,JMAX  Tensor cross sections are computed for all values of the molecular rotational quantum number (lower case $j_i$ and $j_f$ in Eqs. 30 and 31 of the above article) ranging from JMIN to JMAX.

To determine the close-coupled S-matrix elements required by the command TRNPRT, you must carry out a prior calculation with

- The flag [WRSMAT](WRSMAT)= .true.
- JMIN ... JMAX included in the array [JOUT](JOUT)
- IN1 and IN2 included in the array [INDOUT](INDOUT)

The effective cross sections are written in the terminal output and in the file {JOB}{IENERG}.trn.