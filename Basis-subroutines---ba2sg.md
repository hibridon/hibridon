### Collisions between a  <sup>2</sup>&Sigma; Molecule and a Spherical Atom


> BASISTYPE = 2
>
> System subroutine:  sy2sg
>
> Basis subroutine:  hiba02_2sg.F90
>
> Ref. M. H. Alexander, J. Chem. Phys. **76**, 3637 (1982).


#### The definition of the integer system dependent parameters is as follows:

- `NTERM`: (the number of potential surfaces involved:  this should be 2 in this case).  This parameter can not be changed
- `NPAR`: number of lambda-doublet components included
  - if `NPAR` = 1, just &epsilon; = +1 states are included
  - if `NPAR` &#8800; 1, both &epsilon; = +1 and &epsilon; = -1 states are included
- `NRMAX`:  maximum Hund's case (b) rotational quantum number *N*.
- `ISYM`: reflection symmetry of electronic state
  - `ISYM` = +1 for &Sigma;<sup>+</sup> states
  - `ISYM` = -1 for &Sigma;<sup>-</sup> states
- `IGU`: inversion symmetry of states included in expansion
  - `IGU` = +1 for *gerade* states (for heteronuclear molecules `IGU`should be +1)
  - `IGU` = -1 for *ungerade* states
- `ISA`: nuclear permutation symmetry (*s*/*a* label) of states included in the expansion
  If the flag `IHOMO` is .TRUE., then only *s* states will be included if `ISA` = +1 and only *a* states if `ISA` = -1
#### The definition of the real system dependent parameters is as follows:

- `BROT`: rotational constant
- `GSR`: spin-rotation constant
- `DROT` : centrifugal distortion constant
- `HROT`: second centrifugal distortion constant
