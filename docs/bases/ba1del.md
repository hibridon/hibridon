### Collisions between a  <sup>1</sup>&Delta; Molecule and a Spherical Atom


> BASISTYPE = 11
>
> System subroutine:  sy1del
>
> Basis subroutine:  hiba11_1del.F90
>
> Ref. D. G. Sauder, D. Patel-Misra, and P. J. Dagdigian, J. Chem. Phys. **91**, 5316 (1989).


#### The definition of the integer system dependent parameters is as follows:

- `NTERM`: (the number of potential surfaces involved:  this should be 2 in this case).  This parameter can not be changed
- `NPAR`: number of lambda-doublet components included
  - if `NPAR` = 1, just &epsilon; = +1 states are included
  - if `NPAR` &#8800; 1, both &epsilon; = +1 and &epsilon; = -1 states are included
- `JMAX`:  maximum molecular rotational level included in channel basis
- `IGU`: inversion symmetry of states included in expansion
  - `IGU` = +1 for *gerade* states (for heteronuclear molecules `IGU`should be +1)
  - `IGU` = -1 for *ungerade* states
- `ISA`: nuclear permutation symmetry (*s*/*a* label) of states included in the expansion
  If the flag `IHOMO` is .TRUE., then only *s* states will be included if `ISA` = +1 and only *a* states if `ISA ` = -1
#### The definition of the real system dependent parameters is as follows:

- `BROT`: rotational constant
- `Q`: lambda-doubling constat

#### The system dependent input paramaters are as follows:

- line **1**:  `JMAX`, `IGU`, `ISA`, `NPAR`
- line **2**:  `BROT`, `Q`