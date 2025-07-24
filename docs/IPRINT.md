IPRINT is a general print control parameter which can vary from -1 to 2. In some of the [basis](basis) subroutines the amount of diagnostic printing is controlled by IPRINT with the understanding that
  - IPRINT = - 1 : no print
  - IPRINT = 0 : minimal print
  - IPRINT = 1 : moderate diagnostic print
  - IPRINT = 2 : full scale diagnostics

If the flag NOPRIN is set equal to .true., then IPRINT is automatically set equal to -1

At present IPRINT effects printing as follows:

  - BASIS routines: with the flag [BASTST](BASTST) = .TRUE. and IPRINT = 1, only the number of non-zero coupling matrix elements will be printed. To print out the matrix elements themselves, set IPRINT = 2
  - Printing of integral cross sections: For [CS](Coupled-States) calculations, if the flag [PRXSEC](PRXSEC-and-WRXSEC) is set .true. integral cross sections will be printed out only for $J_{\rm initial} and/or $J_{\rm final}$ .le. [NUMAX](NUMAX). If IPRINT = 2, then all the cross sections, even those which are not converged with respect to the sum over the CS projection index $\nu$, are printed.