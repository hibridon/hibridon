### Collision of a <sup>3</sup>&Sigma; diatomic molecule with a closed-shell (<sup>1</sup>&Sigma;) diatomic molecule


> BASISTYPE = 28
>
> System subroutine:  sy3sg1sg
>
> Basis subroutine:  hiba28_3sg1sg.F90
>
> Refs. P. J. Dagdigian, J. Chem. Phys. **150**, 084308 (2019); F. Lique, J. Chem. Phys. **132**, 044311 (2010).


#### The definition of the system dependent parameters is as follows:
- J1MAX: maximum rotational angular momentum of molecule 1
- J2MIN: minimum rotational angular momentum of molecule 2
- J2MAX: maximum rotational angular momentum of molecule 2
- IPOTSY2: cylindrical symmetry of potential. Set to 2 for homonuclear molecule 2, set to 1 for homonuclear molecule 2.

#### The definition of the real system dependent parameters is as follows:

- B1ROT: rotational constant of molecule 1
- D1ROT: centrifugal distortion constant for molecule 1
- GAMMA: spin-rotation constant for molecule 1
- FLMBDA: spin-spin constant for molcule 1
- B2ROT: rotational constant for molecule 2

By setting [`BASTST`](../param/bastst.md)=.TRUE., you can output the values of quantum numbers of the two molecules and the internal energies.