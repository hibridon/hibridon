### Collision of a symmetric top with a linear molecule in a <sup>1</sup>&Sigma; state, using coupled representation of the PES


> BASISTYPE = 21
>
> System subroutine:  systp1sg
>
> Basis subroutine:  hiba21_stp1sg.F90
>
> Ref. C. Rist, M. H. Alexander, and P. Valiron, J. Chem. Phys. **98**, 4662 (1993).


The potential is to be fitted in the following form

$V(R,\theta_1,\phi_1, \theta_2, \phi_2) = \sum_{l_1l_2l,\mu_1\ge0} B_{l_1l_2l\mu_1}(R) (2-\delta_{\mu_1 0}) \sum_{m\ge0} (1+\delta_{m0})^{-1} < l_1m,l2,-m | l0> d^{l_2}_{-m,0}(\theta_2)$

$\times [\ \cos( \mu_1 \phi_1 + m \phi_2 ) d^{l_1}_{\mu_1m}(\theta_1)$

 $+ (-1)^{l_1+l_2+l+\mu_1} \cos( \mu_1 \phi_1 - m \phi_2 ) d^{l_1}_{-\mu_1,m}(\theta_1)\ ]$

The angles are those defined in the above literature. Basis type 9 uses a simpler form of angular expansion, but that expansion is very inefficient in scattering calculations.

Unlike earlier basis routines, this basis routine allows the user to define which terms to include in the angular expansion without limitations other than the total number of terms (which is fixed in Hibridon). This information must be written to the lms array in the mod_bastp1sg module by the loapot subroutine. See pot_stp1sg_qma.f for an example.

Coupled states calculations are NOT supported in this basis routine.

#### The definition of the integer system dependent parameters is as follows:

- IPOTSY: cylindrical symmetry of potential. Presently only symmetric tops with three-fold symmetry are supported. For these system, IPOTSY = 3.

- IOP: bitwise flag for rotational channel basis of different symmetry. The symmetry groups are:
  - Group 0: $k$ is a multiple of 3 and even $j$ for $k=0$, (+) inversion symmetry
  - Group 1: $k$ is a multiple of 3 and odd $j$ for $k=0$, (-) inversion symmetry
  - Group 2: $k$ is not a multiple of 3, (+) inversion symmetry
  - Group 3: $k$ is not a multiple of 3, (-) inversion symmetry
  - Group 4: even $j$ for $k=0$, (-) inversion symmetry
  - Group 5: odd $j$ for $k=0$, (+) inversion symmetry

If group i is to be included in the channel basis, 2i should be added to IOP. Hence, the nuclear permutation symmetry levels of some C3v and D3h molecules are as follows:
  - CH3 (ground vibrational level): ortho-levels have IOP = 2, and para-levels have IOP = 4 or IOP = 8 (the cross sections will be identical using either value).
  - CD3 (ground vibrational level): $A_1$ levels have IOP = 1, $A_2$ levels have IOP = 2, and E-levels have IOP = 4 or IOP = 8 (the cross sections will be identical using either value).
  - NH3: ortho-levels have IOP = 3 (1 + 2, for groups 0 and 1), and para-levels have IOP = 12 (4 + 8, for groups 2 and 3).
  - ND3 (ground vibrational level): $A_1$ levels have IOP = 17, $A_2$ levels have IOP = 34, and E-levels have IOP = 4 or IOP = 8 (the cross sections will be identical using either value).

- J1MAX: the maximum rotational angular momentum in the channel expansion for the symmetric top.

- E1MAX: the maximum energy of a state to be included in the rotational state basis for the symmetric top.

- J2MIN: the minimum rotational angular momentum for the linear molecule.

- J2MAX: the maximum rotational angular momentum for the linear molecule.

- IPOTSY2: step used in iterating the rotational angular momentum for the linear molecule. For H<sub>2&/sub>, IPOTSY2 can be set to 2 so that either para (even J2MIN and J2MAX) or ortho (odd J2MIN and J2MAX) levels are included.

#### The definition of the real system dependent parameters is as follows:

- BROT, CROT: rotational constants of symmetric top

- DROT: rotational constant of diatomic molecule
  ```
   E(jkj2) = BROT j (j + 1) + (CROT - BROT) k<sup>2</sup> + DROT j<sub>2</sub> (j<sub>2</sub> + 1)
  ```
- DELTA: inversion splitting (assumed to be the same for all rotational levels). This energy is added to - inversion levels (groups 1 and 3).

⚠️ At the present time, elastic cross sections are incorrectly calculated using the INTCRS command for channels where both J and J2 are nonzero. Cross sections are correctly calculated using PRINTC.