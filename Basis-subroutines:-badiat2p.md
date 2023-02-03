### Collisions between a  <sup>1</sup>&Delta; Molecule and a Spherical Atom


> BASISTYPE = 15
>
> System subroutine:  sydiat2p
>
> Basis subroutine:  hiba15_diat2p.F90
>
> Ref. M.-L. Dubernet and J. M. Hutson, J. Chem. Phys. **101**, 1939 (1994).


#### The definition of the integer system dependent parameters is as follows:

- NTERM (the number of potential surfaces involved: this should be 3 in this case). This parameter can not be changed

- IOP: *ortho*/*para* label for rotational states of diatomic. If IHOMO =.true. then
  - if IOP = 1: only *para* states included in channel expansion
  - if IOP = -1: only *ortho* states included in channel expansion

- JMAX: maximum rotational quantum number of diatomic.


#### The definition of the real system dependent parameters is as follows:

- BROT: rotational constantof diatomic molecule
- ASO: spin-orbit constant of <sup>2</sup>P atom

