### Collisions between a  <sup>1</sup>&Sigma; Molecule and a Spherical Atom


> BASISTYPE = 1
>
> System subroutine:  sy1sg
>
> Basis subroutine:  hiba01_1sg.F90
>
> Ref.  W. A. Lester, Jr., Meth. Comput. Phys.  **10** , 211 (1971).


#### The definition of the integer system dependent parameters is as follows:

- `NTERM`: (the number of potential surfaces involved:  this should be 1 in this case).  This parameter can not be changed.
- `VIBMN`: vibrational quantum number of lowest included vibrational level
- `VIBMX`: vibrational quantum number of highest included vibrational level\
  ⚠️ to treat rotationally inelastic scattering of a rigid rotor, set VIBMN and VIBMX equal to 1
- `JMIN`:  minimum molecular rotational level included in channel basis
- `JMAX`:  maximum molecular rotational level included in channel basis

#### The definition of the real system dependent parameters is as follows, with values repeated for each vibrational level included:

- `BROT`:    rotational constant
- `DROT`:    centrifugal distortion constant
- `HROT`:    second centrifugal distortion constant\
  ℹ️ E(j) = `BROT` j (j  + 1) -  `DROT`  [ j (j  + 1)] 2  +  `HROT`  [ j (j  + 1)] 3
- `EVIB`: vibrational term (energy of  j  = 0 level) for this vibrational manifold


#### The system dependent input paramaters are as follows:

- line **1**:  `NVIB`, `VIBMN`, `VIBMX`\
  where `NVIB` is the number of different vibrational levels included in the channel basis, and `VIBMN` and `VIBMX` are defined above.\
  ℹ️ Obviously, `NVIB` = `VIBMX` - `VIBMN` + 1
- lines **2**, **5**, **8**,  ...:  `JMIN`, `JMAX`\
  minimum and maximum rotational quantum numbers for  each vibrational manifold
- lines **3**, **6**, **9**,  ...:  `BROT`, `DROT`, `HROT`\
  rotational constants (wavenumbers) for  each vibrational manifold
- lines **4**, **7**, **10**, ...:  `EVIB`\
  vibrational term (wavenumbers) for  each vibrational manifold
