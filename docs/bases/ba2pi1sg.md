### Collisions between a <sup>2</sup>&Pi; molecule and a <sup>1</sup>&Sigma; molecule


> BASISTYPE = 20
>
> System subroutine:  sy2pi1sg
>
> Basis subroutine:  hiba20_2pi1sg.F90
>
> Ref. S. M. Miller and D. C. Clary, J. Chem. Phys. **98**, 1843 (1993); A. R. Offer and M. C. van Hemert, J. Chem. Phys. **99**, 3836 (1993); P. E. S. Wormer et al. J. Chem. Phys. **122**, 244325 (2005).

#### The definition of the integer system dependent parameters is as follows:

- J1MAX: the maximum rotational quantum number, minus 1/2, in the channel expansion for the <sup>2</sup>&Pi; molecule.

- J1MAX is an integer variable. Because we are dealing here with a molecule with half-integer spin, the allowed channel rotational angular momenta run from 1/2 to JMAX+1/2.

- NPAR: number of lambda-doublet components included.
  - If NPAR = 1, just ε = +1 states are included.
  - If NPAR .ne. 1, both ε = +1 and ε = - 1 states are included.

- J2MIN: the minimum rotational angular momentum for the <sup>1</sup>&Sigma; molecule.

- J2MAX: the maximum rotational angular momentum for the <sup>1</sup>&Sigma; molecule.

- IPOTSY2: step used in iterating the rotational angular momentum for the 1Σ molecule.
  - For H<sub>2</sub>, IPOTSY2 can be set to 2 so that either para (even J2MIN and J2MAX) or ortho (odd J2MIN and J2MAX) levels are included.

#### The definition of the real system dependent parameters is as follows:

- BROT: rotational constant of the <sup>2</sup>&Pi molecule

- ASO: spin-orbit constant of the <sup>2</sup>&Pi molecule

- P, Q: lambda-doubling constants of the <sup>2</sup>&Pi molecule

- DROT: rotational constant of the <sup>1</sup>&Sigma; molecule