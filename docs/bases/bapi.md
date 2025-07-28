### Collisions between a <sup>1</sup>, <sup>2</sup>, <sup>3</sup>&Pi; molecule and a spherical Atom


> BASISTYPE = 5
>
> System subroutine:  sypi
>
> Basis subroutine:  hiba05_pi.F90
>
> Ref. M. H. Alexander, J. Chem. Phys. **76**, 5974 (1982); M. H. Alexander, Chem. Phys. **92**, 337 (1985); D. Lemoine, G. C. Corey, M. H. Alexander, and J. Derouard, ibid. **118**, 357 (1987); G. C. Corey and M. H. Alexander, J. Chem. Phys. **85**, 5652 (1986); M. H. Alexander, P. J. Dagdigian, and D. Lemoine, ibid. **95**, 5036 (1991).

#### The definition of the integer system dependent parameters is as follows:

- `NTERM`: (the number of potential surfaces involved: this should be 3 in this case). This parameter can not be changed
- `JMAX`:  the maximum rotational quantum number of diatomic

  For a <sup>2</sup>&Pi; molecule, `JMAX` is the maximum rotational quantum number minus 1/2; so that the allowed channel rotational angular momenta run from 1/2 to `JMAX`+1/2.

- `IGU`: inversion symmetry of states included in expansion
  - `IGU` = +1 for *gerade* states (for heteronuclear molecules `IGU`should be +1)
  - `IGU` = -1 for *ungerade* states
  
   ⚠️ For heteronuclear molecules (`IHOMO` = .False.), `IGU` is ignored but should be set to +1

- `ISA`: nuclear permutation symmetry (*s*/*a* label) of states included in the expansion
  If the flag `IHOMO` is .TRUE., then only *s* states will be included if `ISA` = +1 and only *a* states if `ISA ` = -1

- `NPAR`: number of lambda-doublet components included
  - if `NPAR` = 1, just &epsilon; = +1 states are included
  - if `NPAR` = 1, both &epsilon; = +1 and &epsilon; = -1 states are included

- `IMULT`: spin-multiplicity of &Pi; state (`IMULT` = 2 S + 1)

- `NMAN`: number of spin-orbit manifolds included

  For a <sup>2</sup>&Pi; state, `NMAN` = 1 selects the nominally omega = 1/2 manifold in all other cases, NMAN is set equal to IMULT
#### The definition of the real system dependent parameters is as follows:

- `BROT`: rotational constant 
- `ASO`: spin-orbit constant
- `O`, `P`, `Q`: lambda-doubling constants.

  ⚠️ Not all of these parameters are used for the different multiplicities, see the subroutine lib/hibapi.F90