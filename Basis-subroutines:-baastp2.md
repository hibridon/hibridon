### Collision of an atom with an asymmetric top


> BASISTYPE = 29
>
> System subroutine:  syastp2
>
> Basis subroutine:  hiba29_astp2.F90
>
> Ref. A. Faure et al., ACS Earth Space Sci. **2**, 964 (2019); S. Green, J. Chem. Phys. **64**, 3463 (1976); **67**, 816 (1979).

This version of a basis subroutine for an asymmetric molecule - atom collision applies to a molecule with no symmetry elements, e.g. a chiral molecule. Here, the body-frame z axis is taken to lie along the a inertial axis of the molecule.

#### The definition of the system dependent parameters is as follows:

- NTERM: the number of potential surfaces involved. This parameter can not be changed
- NUMPOT: an index representing the particular potential used. This variable is passed to the POT subroutine
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

By setting [BASTST](BASTST)=.TRUE., you can output the values of $j$, $s$, the values of the prolate and oblate projection quantum numbers [ $k_p$ and $k_o$], and the internal energies, as well as the expansion coefficients.