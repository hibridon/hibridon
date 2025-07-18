### Collisions between of a symmetric top with a linear molecule in a <sup>1</sup>&Sigma; state


> BASISTYPE = 9
>
> System subroutine:  systpln
>
> Basis subroutine:  hiba09_stpln.F90
>
> Ref. C. Rist, M. H. Alexander, and P. Valiron, J. Chem. Phys. **98**, 4662 (1993).

#### The definition of the integer system dependent parameters is as follows:

- `NTERM`: (the number of potential surfaces involved: this should be 2 in this case). This parameter can not be changed
- `NUMPOT`: an index representing the particular potential used. This variable is passed to the POT subroutine
- `IPOTSY`: cylindrical symmetry of potential. If the flag `IHOMO` = .True., only terms with *mu*<sub>1</sub> equal to an integer multiple of `IPOTSY` can be included in the potential. Example: for NH3, `IPOTSY` = 3
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
- `IPOTSY2`: symmetry of potential. if the linear molecule is homonuclear then `IPOTSY` = 2 and only terms with lambda2 even can be included in the potential, else `IPOTSY` = 1
- `J2MAX`: the maximum rotational angular momentum for the linear molecule
- `J2MIN`: the minimum rotational angular momentum for the linear molecule

#### The definition of the real system dependent parameters is as follows:

- `BROT`, `CROT`: rotational constants of symmetric top
- `DROT`: rotational constant of diatomic molecule
*E*(*jkj2*) = `BROT` *j* (*j* + 1) + (`CROT` - `BROT`) *k*<sup>2</sup> + `DROT` *j*<sub>2</sub> (*j*<sub>2</sub> + 1) 
