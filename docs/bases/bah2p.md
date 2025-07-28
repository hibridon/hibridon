### Collisions between an Atom in a <sup>2</sup>P electronic state and a Homonuclear Diatomic


> BASISTYPE = 12
>
> System subroutine:  syh2p
>
> Basis subroutine:  hiba12_h2p.F90
>
> Ref. M.-L. Dubernet and J. M. Hutson, J. Chem. Phys. **101**, 1939 (1994).

#### The definition of the integer system dependent parameters is as follows:

- `NTERM`: (the number of potential surfaces involved: this should be 3 in this case). This parameter can not be changed
- `IOP`: *ortho/para* label for rotational states of diatomic. If `IHOMO` = .True. then
  - if `IOP` = 1: only *para* states are included in channel expansion
  - if `IOP` = -1: only *ortho* states are included in channel expansion
- `JMAX`:  the maximum rotational quantum number of diatomic
#### The definition of the real system dependent parameters is as follows:

- `BROT`: rotational constant of diatomic molecule
- `ASO`: spin-orbit constant of <sup>2</sup>P atom