### Collisions between of an atom with a symmetric top, with no inversion doubling


> BASISTYPE = 18
>
> System subroutine:  systp1
>
> Basis subroutine:  hiba18_stp1.F90
>
> Ref. S. Green, J. Chem. Phys. **64**, 3463 (1976); M. H. Alexander, J. Chem. Phys. **77**, 1855 (1982).

#### The definition of the integer system dependent parameters is as follows:

- `NTERM`: (the number of potential surfaces involved: this should be 1 in this case). This parameter can not be changed
- `NUMPOT`: an index representing the particular potential used. This variable is passed to the POT subroutine
- `IPOTSY`: cylindrical symmetry of potential. If the flag `IHOMO` = .True., only terms with *mu* equal to an integer multiple of `IPOTSY` can be included in the potential. Example: for NH3, `IPOTSY` = 3
- `IOP`: *ortho/para* label for rotational states of diatomic. If `IHOMO` = .True. then
  - if `IOP` = 1: only *para* states are included in channel expansion
  - if `IOP` = -1: only *ortho* states are included in channel expansion
- `NINV`:  number of inversion doublets included:
  - if `NINV` = +1, only + inversion levels included
  - if `NINV` = -1, only - inversion levels included
  - if `NINV` = 2, both inversion levels included

- `JMAX`: the maximum rotational angular momentum included in the channel expansion
- `EMAX`: the maximum energy of a state to be included in the rotational state basis

#### The definition of the real system dependent parameters is as follows:

- `BROT`, `CROT`: rotational constants of symmetric top
*E*(*jk) = `BROT` *j* (*j* + 1) + (`CROT` - `BROT`) *k*<sup>2</sup>
