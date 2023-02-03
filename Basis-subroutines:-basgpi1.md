### Collisional mixing of <sup>2</sup>&Pi; and <sup>2</sup>&Sigma; states of a diatomic molecule induced by a spherical atom (no spectroscopic perturbations)


> BASISTYPE = 19
>
> System subroutine:  sysgpi1
>
> Basis subroutine:  hiba19_sgpi1.F90
>
> Ref. M. H. Alexander and G. C. Corey, J. Chem. Phys. **84**, 100 (1986); H.-J. Werner, B. Follmeg, M. H. Alexander, and D. Lemoine, ibid. **91**, 5425 (1989).


The vibrational levels included involve one level in the $^2\Sigma$ state and one or more levels in the $^2\Pi$ state.

#### The definition of the integer-system dependent parameters is as follows:

- NTERM (the number of potential surfaces involved: this can be 4, 7, 10. etc. here and must be consistent with the POT subroutine). This parameter cannot be changed
- ISYM: if ISYM =+1, then the electronic symmetry of the $^2\Sigma$ state is sigma-plus if ISYM = -1, then the electronic symmetry is sigma-minus
- ISA: s/a symmetry index, if the molecule is homonuclear (IHOMO=T) then, if ISA=+1 then only the s-levels (both for sigma and pi) are included in the basis, if ISA=-1, then only the a-levels are included
- IGUSG: permutation inversion symmetry of $^2\Sigma$ electronic state (for homonuclear molecules, IHOMO=T)
  - IGUSG = 1 for gerade states
  - IGUSG = -1 for ungerade states
- NMAXSG: highest rotational level for the 2 state
- NPARSG: number of spin doublets in $^2\Sigma$ state included
  - NPARSG = 2 will ensure both spin doublets
  - NPARSG = 1, just  = 1 doublets
  - NPARSG = -1, just  = -1 doublets
- IGUPI: permutation inversion symmetry of $^2\Pi$ electronic state (for homonuclear molecules, IHOMO=T)
  - IGU = 1 for gerade states
  - IGU = -1 for ungerade states
- NPARPI: number of $^2\Pi$ symmetry doublets included
  - NPARPI = 2 will ensure both lambda doublets
  - NPARPI = 1, just  = 1 doublets
  - NPARPI = -1, just  = -1 doublets
- NUMVPI: number of $^2\Pi$ vibrational levels to be included

- IVPI: array of vibrational quantum number for $^2\Pi$ vibrational levels. Must be consistent with the potential routine

- JMAX: array of maximum rotational angular momenta for each $^2\Pi$ channel

In each spin-orbit manifold with convention omega .le. j .le. jmax+0.5

 ⚠️ for heteronuclear molecules both IGUPI and IGUSG should be +1