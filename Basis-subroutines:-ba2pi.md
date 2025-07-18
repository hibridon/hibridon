### Collisions between a  <sup>2</sup>&Pi; Molecule and a Spherical Atom


> BASISTYPE = 3
>
> System subroutine:  sy2pi
>
> Basis subroutine:  hiba03_2pi.F90
>
> Ref. M. H. Alexander, J. Chem. Phys. **76**, 5974 (1982); M. H. Alexander, Chem. Phys. **92**, 337 (1985); G. C. Corey and M. H. Alexander, J. Chem. Phys. **85**, 5652 (1986).


#### The definition of the integer system dependent parameters is as follows:

- `NTERM`: (the number of potential surfaces involved:  this should be 2 in this case).  This parameter can not be changed
- `NPAR`: number of lambda-doublet components included
  - if `NPAR` = 1, just &epsilon; = +1 states are included
  - if `NPAR` &#8800; 1, both &epsilon; = +1 and &epsilon; = -1 states are included
- `JMAX`:  maximum rotational quantum number **minus 1/2**
- `IGU`: inversion symmetry of states included in expansion
  - `IGU` = +1 for *gerade* states (for heteronuclear molecules `IGU`should be +1)
  - `IGU` = -1 for *ungerade* states
- `ISA`: nuclear permutation symmetry (*s*/*a* label) of states included in the expansion
  If the flag `IHOMO` is .TRUE., then only *s* states will be included if `ISA` = +1 and only *a* states if `ISA ` = -1
#### The definition of the real system dependent parameters is as follows:

- `BROT`: rotational constant
- `ASO` : spin-orbit constant
- `P`, `Q`: lambda-doubling constats
