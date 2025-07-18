The OPTIMIZE command allows the user to determine easily optimal values for the following system independent [parameters](Parameters) which govern the accuracy of the Hibridon code:

[FSTFAC](FSTFAC), [RCUT](RCUT), [RENDAI](RSTART,-RENDLD,-and-RENDAI), [RENDLD](RSTART,-RENDLD,-and-RENDAI), [RINCR](TOLAI-and-RINCR), [RSTART](RSTART,-RENDLD,-and-RENDAI), [SPAC](SPAC), [TOLAI](TOLAI-and-RINCR)

The command line syntax is
```
OPTIMIZE,{pname},IVAL,FVAL,FACT,STEP,ACCAV,ACCMX,THRESH
```
(upper or lower case), where
- {pname}:  Name of the parameter from the above list which is to be varied

- IVAL:  Initial value of this parameter

- FVAL:  Final value of this parameter

- FACT,STEP:  Variables which control the step size for the optimization.

The value of the parameter PNAME is updated at each step following the equation
```
new value = (old value) * FACT + STEP
```

- ACCAV, ACCMAX, THRESH:  The optimization is stopped when, for two subsequent runs, both the average percentage difference between the S matrices (or the moduli of the Smatrices) is smaller than ACCAV and the maximum percentage difference between the S matrices (or between the moduli of the Smatrices) is smaller than ACCMAX. In this comparison all S-matrix elements of magnitude less than THRESH are neglected. If THRESH is positive, or zero, then the comparison involves the moduli of the S-matrix elements. If THRESH is negative, then the comparison involves both the real and imaginary parts of the S matrix.

The optimization is performed only for $J_{\rm tot}$ = [JTOT1](JTOT)

Values for {pname}, IVAL, and FVAL must be given. If values for the other variables are omitted (see below), then default values are assigned as follows:

- ACCAV = 1
- ACCMAX = 5
- THRESH = 1.E-5

If FACT and STEP are not defined, then FACT is set equal to
- 2, if FVAL > IVAL
- 0.5, if FVAL < IVAL

 Values of the various input parameters may be omitted, but their place must be preserved by commas if any of the subsequent parameters are to be given values other than their default values.