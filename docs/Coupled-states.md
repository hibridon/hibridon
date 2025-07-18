## Introduction

In the Coupled States approximation [1-5] the diagonal centrifugal matrix of the centrifugal barrier in the [close-coupled](Close-coupled-equations) formulation

$(l^2)_{ii} = \frac{h^2}{2mR^2} l_i (l_i + 1)$, (1)

is replaced by a constant diagonal matrix

$(^l2)_{ii} = \frac{h^2}{2mR^2}  \bar{L}(\bar{L}+ 1)$, (2)

where $\bar{L}$ (designated L-bar) is some average orbital angular momentum, then the CC equations become block-diagonal in the index  which is the body-frame projection of $J$. Equivalently, the block diagonalization can be achieved by neglecting all Coriolis coupling in a body-frame expansion of the total scattering wavefunction, and, in addition, approximating the diagonal terms in the body-frame expansion of the orbital angular momentum

$\mathbf{l} = \mathbf{J}_{\rm tot} - \mathbf{J}$

by Eq. (2).

Because the computer time required for solution of the CC equations scales as $N_{\rm ch}^3$, this block diagonalization leads to a substantial reduction in the total time required in a scattering calculation. Moreover, for problems in which the potential is not long-ranged, many tests have shown the the CS approximation is remarkably accurate.

## Cross Sections

Within the CS formulation, the integral cross section for a transition from an initial state $i$ to a final state $f$ is given by
 
$\sigma_{if} = \frac{\pi}{(2j_i+1)k_i^2} \sum_{\bar{L}, \nu} |\delta_{if} - S_{if}^{\bar{L}, \nu}|^2 $(3)

where $k_i$ is the wavevector of the initial state and the sum runs over:
- all values of the effective orbital angular momentum 
- all values of the projection quantum number  such that | $\nu$ | .le. min { $\mathbf{J}_i$, $\mathbf{J}_f$ }

Note that in the CS formulation, the S matrix is indexed in both $\bar{L}$ and $\nu$. Since the S matrix is independent of the sign of $\nu$ , the summation in Eq. (3) can be further restricted to positive, definite values of $\nu$ 
At present the Hibridon code contains no provision for the determination of differential cross sections with the CS approximation.


***

References

1. P. McGuire and D. J. Kouri, J. Chem. Phys. **60**, 2488 (1974).
2. P. McGuire, Chem. Phys. **13**, 81 (1976).
3. R. T. Pack, J. Chem. Phys. **60**, 633 (1974).
4. D. J. Kouri in Atom-Molecule Collision Theory: A Guide for the Experimentalist, edited by R. B. Bernstein (Plenum, New York, l979) p. 301. 
5. A. E. DePristo and M. H. Alexander, Chem. Phys. 19, **181** (1977).