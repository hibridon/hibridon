### Collisions of two unlike a closed-shell (<sup>1</sup>&Sigma;) diatomic molecule


> BASISTYPE = 25
>
> System subroutine:  sy1sg1sg
>
> Basis subroutine:  hiba25_1sg1sg.F90
>
> Ref. S. Green, J. Chem. Phys. **62**, 2271 (1975).

#### The definition of the system dependent parameters is as follows:

- J1MAX: maximum rotational angular momentum of molecule 1

- J2MIN: minimum rotational angular momentum of molecule 2

- J2MAX: maximum rotational angular momentum of molecule 2

- IPOTSY: cylindrical symmetry of potential. Set to 2 for homonuclear molecule 2, set to 1 for homonuclear molecule 2.

#### The definition of the real system dependent parameters is as follows:

- B1ROT: rotational constant of molecule 1

- D1ROT: centrifugal distortion constant for molecule 1

- B2ROT: rotational constant for molecule 2

By setting [BASTST](BASTST)=.TRUE., you can output the values of quantum numbers of the two molecules and the internal energies.