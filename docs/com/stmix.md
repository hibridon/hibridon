# STMIX

The command STMIX allows you to calculate cross sections for states with mixed singlet-triplet electronic character, as in the case of mixed CH2(3B1,1A1) levels. It is assumed that the perturbation connects a single level in the singlet manifolc with a single level in the triplet manifold. The command line is
```
STMIX,JOB1,IENERG1,JOB2,IENERG2,DELE,EMAX,ISTATS,ISTATT,WSO
```
where the input parameters are:
- JOB1,IENERG1  These parameters locate the S-matrix for rotational levels in the singlet manifold. This matrix should be computed using the basis subroutine BAASTP. The program searches for an S-matrix in the file Job1{ienerg1}.smt, where IENERG1 is an integer (IENERG1 = 1, 2, ... etc.)

- JOB2,IENERG2  These parameters locate the S-matrix for rotational levels in the triplet manifold. This matrix should be computed in a spin-free calculation using the basis subroutine BACH2X. The program searches for an S-matrix in the file Job2{ienerg1}.smt, where IENERG2 is an integer (IENERG2 = 1, 2, ... etc.)

- DELE   The difference in energy between the origins of the singlet and triplet manifold, computed as ETOT(singlet) − ETOT(triplet). The subroutine checks that EDELE equals this energy difference, to within 0.1 cm−1.

- EMAX  The maximum energy of a state to be included in the rotational state basis.

- ISTATS  The state number of the perturbed level in the singlet manifold. The state number can be obtained by running STMIX with WSO set to zero.

- ISTATT  The state number of the perturbed level in the triplet manifold. The state number can be obtained by running STMIX with WSO set to zero. It should be noted that the triplet state list includes the spin states with each value of the rotational angular momentum n. The subroutine checks that the total angular momentum j is the same for the perturbed singlet and triplet levels.

- WSO  The value of the spin-orbit matrix element between the two perturbed levels. As noted above, if WSO is set to zero, the channel list will be output, but cross sections will not be computed.