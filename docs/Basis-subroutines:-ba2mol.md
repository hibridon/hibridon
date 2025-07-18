### Collisions between of two <sup>1</sup>&Sigma; Diatomics


> BASISTYPE = 8
>
> System subroutine:  sy2mol
>
> Basis subroutine:  hiba8_2mol.F90
>
> Ref. S. Green, J. Chem. Phys. **62**, 2271 (1975); A. E. DePristo and M. H. Alexander, ibid. **66**, 1334 (1977); M. H. Alexander, ibid. **73**, 5135 (1980)


#### The definition of the integer system dependent parameters is as follows:

- `NTERM`: (the number of potential surfaces involved: this should be 1 in this case).  This parameter can not be changed
- `NSYM`:  interchange symmetry of the included channels
- `NUMJ`:  number of different `J1`/`J2` pairs in channel list
- `J1`/`J2`:  vector of length `NUMJ` containing packed index of rotational quantum number of each molecule with *j*<sub>1</sub> in the tens digit and j<sub>2</sub> in the ones digit.

For example, if `NUMJ` = 3, the program expects three `J1`/`J2` values. Suppose these values are, successively, 00, 11, and 02. The calculation will then include all channels with the following pairs of molecular rotational angular momenta:
  - *j*<sub>1</sub> = 0, *j*<sub>2</sub> = 0
  - *j*<sub>1</sub> = 1, *j*<sub>2</sub> = 1
  - *j*<sub>1</sub> = 0, *j*<sub>2</sub> = 2

Note that because of interchange symmetry, the last entry will automatically include all channels with *j*<sub>1</sub> = 2, *j*<sub>2</sub> = 0
#### The definition of the real system dependent parameters is as follows:

- `BROT`: rotational constant
- `DROT`: centrifugal distortion constant
- `HROT`: second centrifugal distortion constant

  *E*(*j*) = *B* *j* (*j* + 1) - *D* [*j* (*j* + 1)]<sup>2</sup> + *H* [*j* (*j* + 1)]<sup>3</sup>

