The command TENXSC allows you to calculate the tensor cross sections defined by Eqs. (30) and (31) of B. Follmeg, P. Rosmus, and H.-J. Werner, J. Chem. Phys. **93**, 4687 (1990). The command line is
```
TENXSC,JOB,IENERG,IFRAME,LAMMAX,KMAX,IN1,IN2,JTOTMX,JMIN,JMAX
```
where the input parameters are:
- JOB,IENERG  The program searches for an S-matrix in the file Job{ienerg}.smt, where IENERG is an integer (IENERG = 1, 2, ... etc.)
The default values of JOB is JOB and of INERG is 1
- IFRAME   This parameter defines the quantization frame employed in the tensor cross section calculations (angle-integrated cross sections, except where indicated otherwise below).
  - IFRAME=0: The quantization axis is defined as the laboratory Z-axis, with an isotropic relative velocity distribution.
  - IFRAME=1: The quantization axis lies along the initial relative velocity vector (collision frame).
  - IFRAME=2: The initial level has quantization axis along the initial relative velocity vector, and the final level has quantization axis along the final velocity vector (helicity frame). In this case, m-resolved differential cross sections are computed over a grid of scattering angles (0.5 deg spacing) for the (J1MIN,IN1) -> J1MIN,IN1) elastic transition, and then differential tensor cross sections are computed. The differential tensor cross sections are numerically integrated to determine integral tensor cross sections.

- LAMMAX   The maximum value of the anisotropic term in the velocity distribution (lower case lambda in the above article). Tensor cross sections for all even values of lambda from 0 to LAMMAX are calculated.  ⚠️ At present only LAMMAX=0 is operational!

- KMAX   The maximum value of the multipole order for which tensor cross sections are calculated. Cross sections for all values of $K_i$ and $K_f$ .le. KMAX are calculated.

- IN1,IN2  The initial and final values of the additional index

- JTOTMX  The maximum value included in the summation over total angular momentum (upper case J in Eqs. 30 and 31 of the above article)

- JMIN,JMAX  Tensor cross sections are computed for all values of the molecular rotational quantum number (lower case $j_i$ and $j_f$ in Eqs. 30 and 31 of the above article) ranging from JMIN to JMAX.

To determine the S-matrix elements required by the command TENXSC, you must carry out a prior calculation with

- the flag [WRSMAT](WRSMAT)= .true.
- JMIN ... JMAX included in the array [JOUT](JOUT)
- IN1 and IN2 included in the array [INDOUT](INDOUT)