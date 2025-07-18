### Collision of an Atom with a CH<sub>2</sub>(X <sup>3</sup>B<sub>1</sub>) (0,v<sub>2</sub>,0) bender vibrational level


> BASISTYPE = 17
>
> System subroutine:  sych2x
>
> Basis subroutine:  hiba17_ch2x.F90
>

#### The definition of the system dependent parameters is as follows:

- `NTERM`: the number of potential surfaces involved. This parameter should equal the number of MU terms in the angular expansion of the potential.. This parameter can not be changed

- `NUMPOT`: an index representing the particular potential used. This variable is passed to the POT subroutine

- `IPOTSY`: cylindrical symmetry of potential. The variables (theta,phi) describing the angular expansion of the potential should be defined with the *a* inertial axis defined as the body-frame *z* axis and, if possible, the *xz* plane as a plane of symmetry of the molecule (in this case, `POTSY` = 2). If the flag `IHOMO` = .True., only terms with `LAMBDA` + `MU` equal to an integer multiple of `IPOTSY` can be included in the potential. Example: for H<sub>2</sub>O, `IPOTSY` = 2

- `IOP`: ortho/para label for molecular states of the asymmetric top. If `IHOMO` =.True. then
  - if `IOP` = 1: only para states included in channel expansion
  - if `IOP` = −1: only ortho states included in channel expansion
- `IVBEND`: the v2 bending vibrational quantum number [0 <= IVBEND <= 3]
- `JMAX`: the maximum rotational angular momentum included in the channel expansion
- `EMAX`: the maximum energy of a state to be included in the rotational state basis


The CH<sub>2</sub>(X <sup>3</sup>B<sub>1</sub>) molecule has a low barrier to linearity, and its rotational energies are not at all well described by the standard rotational energy formulas. Bunker and Jensen have developed a Morse oscillator-rotating bender internal dynamics Hamiltonian (MORBID) to compute the rovibrational energies in this electronic state [see P. Jensen and P. R. Bunker, J. Chem. Phys. **89**, 1327 (1988); P. R. Bunker, P. Jensen, W. P. Kraemer, and R. Beardworth, J. Chem. Phys. **85**, 3724 (1986)]. Rather than compute the rotational energies in a given vibrational level, here we look up the energies in a table provided by Jensen [P. Jensen, private communication (2010)].
The CH<sub>2</sub>(X <sup>3</sup>B<sub>1</sub>) state is well described by Hund's case (b) coupling, in which the electron spin [S=1] is weakly coupled to the molecular frame. As a result, the rotational rotational angular momentum n is coupled to the electron to yield states of total angular momentum *j* = *n* − 1, *n*, *n* + 1 [if *n* >= 1]. This subroutine sets up a calculation of spin-free cross sections, for transitions between states of differing *n*.

Because of the large effective *A* rotational constant in the bender levels, the rotational wave functions are well approximated by symmetric top wave functions. The latter can be expressed in a symmetrized basis as [S. Green, J. Chem. Phys. **64**, 3463 (1976)]

| *n k m s* > = [2(1+&delta;<sub>k0</sub>)]−1/2 [ | *n k m > + s* | *n −k m* >]
where s is the symmetry index (+1 or −1). By setting `BASTST` = .TRUE., you can output the values of *n*, *s*, the values of the asymmetric top prolate and oblate projection quantum numbers [*k*<sub>p</sub> and *k*<sub>o</sub>], and the internal energies.