### Collisions of an asymmetric top with a closed-shell linear molecule


> BASISTYPE = 30
>
> System subroutine:  syastp3
>
> Basis subroutine:  hiba30_astp3.F90
>
> Refs. T. R. Phillips, S. Maluendes, and S. Green, J. Chem. Phys. **102**, 6024 (1995); P. Valiron et al., J. Chem. Phys. **129**, 134306 (2008).

This basis subroutine pertains to collisions of an asymmetric top molecule of $C_{2v}$ or $C_s$ symmetry with a closed-shell diatomic molecule. 
For a $C_{2v}$ molecule, the body-frame z axis is taken to lie along the $C_{2}$ axis,, following Green's convention.

For a molecule with $C_{v}$ symmetry, the body-frame z axis is taken to lie along the an inertial axis.


#### The definition of the integer system dependent parameters is as follows:

- NTERM: the number of potential surfaces involved. This parameter can not be changed
- NUMPOT: an index representing the particular potential used. This variable is passed to the POT subroutine
- IOP: ortho/para label for molecular states of the asymmetric top of $C_{2v}$ symmetry. If [IHOMO](IHOMO) =.true. then
  - if IOP = 1: only para states included in channel expansion
  - if IOP = −1: only ortho states included in channel expansion
  IOP is not needed if [IHOMO](IHOMO) =.false. (for a molecule of $C_{s}$ symmetry, but a number should be entered for this parameter nevertheless.

- JMAX: the maximum rotational angular momentum included in the channel expansion

- IPOTSY2: symmetry of potential. Set equal to 2 for homonuclear molecule 2, or 1 for heteronuclear molecule.

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