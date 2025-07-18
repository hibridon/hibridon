### Collision of a spherical top molecule with an atom 


> BASISTYPE = 24
>
> System subroutine:  sysphtp
>
> Basis subroutine:  hiba25_sphtp.F90
>
> Ref. J. M. Hutson and A. E. Thornley, J. Chem. Phys. **100**, 2505 (1994); T. G. A. Hiejmen, T. Korona, R. Moszynski, P. E. B. Wormer, and A. van der Avoird, J. Chem. Phys. **107**, 902 (1997).

#### The definition of the integer system dependent parameters is as follows:

- NTERM: the number of potential surfaces involved. This parameter should equal the number of MU terms in the angular expansion of the potential. This parameter can not be changed

- NUMPOT: an index representing the particular potential used. This variable is passed to the POT subroutine

- IOP: nuclear spin modification label for the molecular states:
  - if IOP = 1: only A levels included in channel expansion
  - if IOP = 2: only E levels included in channel expansion
  - if IOP = 3: only F levels included in channel expansion

- JMAX: the maximum rotational angular momentum included in the channel expansion

#### The definition of the real system dependent parameter is as follows:

- BROT: rotational constant of the sphericxal top