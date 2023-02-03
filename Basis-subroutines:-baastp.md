### Collisions of an Atom with an Asymmetric Top


> BASISTYPE = 16
>
> System subroutine:  syatp
>
> Basis subroutine:  hiba02_2sg.F90
>
> Ref. B. J. Garrison et al., J. Chem. Phys. *65*, 2193 (1976); S. Green, J. Chem. Phys. *64*, 3463 (1976); 67, 816 (1979).


#### The definition of the integer system dependent parameters is as follows:

- `NTERM`: the number of potential surfaces involved. This parameter should equal the number of `MU` terms in the angular expansion of the potential. This parameter can not be changed
- `NUMPOT`: an index representing the particular potential used. This variable is passed to the `POT` subroutine
- `IPOTSY`: cylindrical symmetry of potential. The variables (&theta;,&phi;) describing the angular expansion of the potential should be defined with the a inertial axis defined as the body-frame *z* axis and, if possible, the *xz* plane as a plane of symmetry of the molecule (in this case, `IPOTSY` = 2). If the flag `IHOMO` = .true., only terms with `LAMBDA` + `MU` equal to an integer multiple of `IPOTSY` can be included in the potential. Example: for H<sub>2</sub>O, `IPOTSY` = 2
- `IOP`: *ortho/para* label for molecular states of the asymmetric top. If `IHOMO` =.true. then
  - if `IOP` = 1: only *para* states are included in channel expansion
  - if `IOP` = -1: only *ortho* states are included in channel expansion
- `JMAX`:  the maximum rotational angular momentum included in the channel expansion
- `EMAX`: the maximum energy of a state to be included in the rotational state basis
#### The definition of the real system dependent parameters is as follows:

- `AROT`, `BROT`, `CROT`: rotational constants of the asymmetric top

---

The rotational eigenfunctions of an asymmetric top are expanded in a symmetrized symmetric top basis as: [S. Green, J. Chem. Phys. **64**, 3463 (1976)]

![image](https://user-images.githubusercontent.com/20931653/171144778-1c6802bc-6d98-4d72-8ad0-b771662b34db.png)

where *s* is the symmetry index (+1 or −1). In this basis, the asymmetric top hamiltonian block diagonalizes into 4 groups: 
  - (1) *k* even, *s* = +1
  - (2) *k* even, *s* = −1
  - (3) *k* odd, *s* = +1
  - (4) *k* odd, *s* = −1

The expansion coefficients are stored in the array c as:
  - *c*(*k* = 0), *c*(*k* = 2), ... *c*(*k* = *j*), for even *k*
  - *c*(*k* = 1), *c*(*k* = 3), ... *c*(*k* = *j*), for odd *k*

By setting `BASTST`=.TRUE., you can output the values of *j*, *s*, the values of the prolate and oblate projection quantum numbers [*k<sub>p</sub>* and *k<sub>o</sub>*], and the internal energies, as well as the expansion coefficients.
