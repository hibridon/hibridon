### Collisions of an atom with an asymmetric top


> BASISTYPE = 27
>
> System subroutine:  syastp1
>
> Basis subroutine:  hiba27_astp1.F90
>
> Ref. B. J. Garrison et al., J. Chem. Phys. **65**, 2193 (1976); S. Green, J. Chem. Phys. **64**, 3463 (1976); **67**, 816 (1979).

Here, the body-frame z axis is taken to lie along the $C_2$ axis of a molecule with $C_{2v}$ symmetry, following Green's convention.

#### The definition of the integer system dependent parameters is as follows:

- NTERM: the number of potential surfaces involved. This parameter can not be changed

- NUMPOT: an index representing the particular potential used. This variable is passed to the POT subroutine

- IPOTSY: cylindrical symmetry of potential. The variables (theta,phi) describing the angular expansion of the potential should be defined with the a inertial axis defined as the body-frame z axis and, if possible, the xz plane as a plane of symmetry of the molecule (in this case, POTSY = 2). If the [flag](../Flags.md) [`IHOMO`](../param/ihomo.md) = .true., only terms with LAMBDA + MU equal to an integer multiple of IPOTSY can be included in the potential. Example: for H2O, IPOTSY = 2

- IOP: ortho/para label for molecular states of the asymmetric top. If [`IHOMO`](../param/ihomo.md) =.true. then
  - if IOP = 1: only para states included in channel expansion
  - if IOP = −1: only ortho states included in channel expansion

- JMAX: the maximum rotational angular momentum included in the channel expansion

#### The definition of the real system dependent parameters is as follows:

- AROT, BROT, CROT: rotational constants of the asymmetric top

- EMAX: the maximum energy of a level to be included in the rotational state basis

***

The rotational eigenfunctions of an asymmetric top are expanded in a symmetrized symmetric top basis as [S. Green, J. Chem. Phys. **64**, 3463 (1976)]

$| jkms > = [2(1+ \delta k_0)]^{1/2} [|jkm> + s |j-km>] $

where $s$ is the symmetry index (+1 or −1). In this basis, the asymmetric top hamiltonian block diagonalizes into 4 groups: (1) $k$ even, $s$ = +1, (2) $k$ even, s = −1, (3) $k$ odd, s = +1, and (4) $k$ odd, s = −1.

The expansion coefficients are stored in the array c as:

- $c(k=0)$, $c(k=2)$, ... $c(k=j)$, for even $k$
- $c(k=1)$, $c(k=3)$, ... $c(k=j)$, for odd $k$

By setting [`BASTST`](../param/bastst.md)=.TRUE., you can output the values of $j$, $s$, the values of the prolate and oblate projection quantum numbers [ $k_p$ and $k_o$], and the internal energies, as well as the expansion coefficients.