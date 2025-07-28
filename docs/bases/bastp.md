### Collisions between of an atom with a symmetric top


> BASISTYPE = 6
>
> System subroutine:  systp
>
> Basis subroutine:  hiba06_stp.F90
>
> Ref. S. Green, J. Chem. Phys. **64**, 3463 (1976); **67**, 816 (1979).

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
- `KMAX`: the maximum projection quantum number included
- `JMAX0`:  the maximum rotational angular momenta or the k = 0 stack
- `JMAX1`:  the maximum rotational angular momenta or the k = 1 stack
- `JMAX2`:  the maximum rotational angular momenta or the k = 2 stack
- `JMAX3`:  the maximum rotational angular momenta or the k = 3 stack
- `JMAX4`:  the maximum rotational angular momenta or the k = 4 stack
- `JMAX5`:  the maximum rotational angular momenta or the k = 5 stack
- `JMAX6`:  the maximum rotational angular momenta or the k = 6 stack
- `JMAX7`:  the maximum rotational angular momenta or the k = 7 stack
- `JMAX8`:  the maximum rotational angular momenta or the k = 8 stack
- `JMAX9`:  the maximum rotational angular momenta or the k = 9 stack
- `JMAX10`:  the maximum rotational angular momenta or the k = 10 stack
- `JMAX11`:  the maximum rotational angular momenta or the k = 11 stack

#### The definition of the real system dependent parameters is as follows:

- `BROT`, `CROT`: rotational constants of symmetric top
*E*(*jk) = `BROT` *j* (*j* + 1) + (`CROT` - `BROT`) *k*<sup>2</sup>
