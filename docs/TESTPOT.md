The command TESTPOT allows you to print out the value of either
- The potential V(R,) for a given value of the arguments
- The expansion coefficients V(R) for a given value of R.

You will be queried which type of output you wish.

⚠️ Meaningful values of $V(R,\theta)$ will be given only if the potential used is expanded as

$V(R,\theta) = V_{00}(R) + \sum_\lambda V_\lambda(R) \times P_\lambda[cos(\theta)]$

where $\lambda$ = 1, 2, ..., $P_\lambda$ is an ordinary Legendre polynomial and the sum extends from LAMMIN(1) to min[LAMMAX(1), 6].
In the case of a homonuclear molecule ([IHOMO](IHOMO-and-TWOMOL)=.true.), only even terms are included in the sum.

The functions $V_{00}(R)$ and $V_\lambda(R)$ are calculated by the subroutine POT.

After entering the command TESTPOT you will be prompted for the value of nlam, which is the maximum number of anisotropic terms in the expansion of the potential.