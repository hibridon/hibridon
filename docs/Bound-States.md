## Introduction
In an entirely equivalent manner to the description of molecular collisions, one can make a [Close-Coupled (CC)](Close-coupled-equations) expansion for wavefunctions of bound states of weakly bound complexes. For more detail we refer the reader to an excellent review by Hutson.[1] The wavefunction of the complex is expanded in a complete set of internal states of the system, usually constructed as direct products of the internal states of one (or both) fragments, multiplied by angular functions which describe the rotation of one collision partner about the other. As in the [`CC`](Close-coupled-equations) description of scattering we designate these internal states as (r), where r designates the internal coordinates. Each internal state is called a channel. The full wavefunction of the complex is written as:

$\Psi(\mathbf{R}, \mathbf{r}) = \sum \mathbf{F}(\mathbf{R}) \times \Phi(\mathbf{r})$    (1)

The nth column of the \mathbf{F}(\mathbf{R}) matrix defines the expansion coefficients for the nth bound state.
Substitution of expansion (1) into the Schrödinger equation, premultiplication by one of the internal states, and integration over r gives rise to a set of coupled ordinary differential equations for the expansion coefficients \mathbf{F}(\mathbf{R}), which are identical to the [`CC`](Close-coupled-equations) description of molecular scattering.

The general structure of these coupled second-order differential equations is expressed by the matrix equation:

$\left[ \mathbf{1} \frac{\partial^2}{\partial R^2} + \mathbf{W}(R) \right] \mathbf{F}(R) = 0$.   (2)

Here $\mathbf{1}$ designates the identity matrix, $R$ is the interparticle distance, and the matrix $\mathbf{W}(R)$ is given by

$\mathbf{W}(R) = \mathbf{k}^2 - \mathbf{l}^2 - \frac{2m}{h^2} \mathbf{V}(R)$,   (3)

where $h$ designates Planck's constant divided by 2, $m$ is the reduced mass of the collision system and $\mathbf{k}^2$ and $\mathbf{l}^2$ designate, respectively, the (diagonal) matrices of the wavevector and the relative orbital angular momentum of the collision partners. We have

$(k^2)_{ii} = \frac{2m}{h^2}(E-e_i)$,    (4)

where $E$ is the total energy and $e_i$ is the internal energy of the ith channel. Also, we have

$(l^2)_{ii} = \frac{h^2}{2mR^2} l_i (l_i + 1)$,   (5)

where l_i is the relative orbital angular momentum in the ith channel. In Eq. (3) the matrix $\mathbf{V}(R)$ is the (full, symmetric) matrix of the coupling potential.

Diagonalization of the $\mathbf{W}(R)$ matrix yields the diagonal matrix of adiabatic wavevectors $\mathbf{k}(R)$. The eigenvectors define the locally adiabatic states, which are transformations of the internal states used to expand the scattering wavefunction, $\Phi(\mathbf{r})$. If $\mathbf{C}(R)$ designates the matrix of eigenvectors, column ordered, then the diagonal matrix of adiabatic energies is defined as

$e(R) = \mathbf{C}(R)^T \mathbf{V}(R) \mathbf{C}(R)$.   (6)

## Energies and Wavefuntions

In the adiabatic bender treatment of bound states,[2, 3] these adiabatic energies define a set of one-dimensional potentials in R. The wavefunctions and energies of these one-dimensional potentials can be obtained by any standard numerical integration technique, such as the Numerov method.[4] These wavefunctions are the adiabatic-bender approximations to the full wavefunctions of the bend-stretch states of the complex.
Alternatively, the energies and wavefunctions can be obtained by a dual numerical propagation of the the matrix of solutions F(R). Propagation is outward from a value of the interparticle distance $R = R_{\rm start}$, which lies well inside the innermost classical turning point, and inward from $R = R_{\rm end}$, which lies well outside the outermost classical turning point. The outwardly and inwardly propagated solution matrices are matched at a distance $R = R_{\rm match}$. Only at the allowed energies of the system will you be able to match, without discontinuity, both the wavefunction and its first derivative.

Great numerical stability is obtained by propagation of the logarithmic derivative of the solution matrix F(R), namely

$\mathbf{Y}(R) = \mathbf{F}'(R)\mathbf{F}(R)- 1$  (7)

rather than the solution matrix itself. The matching condition is [5]

$|\mathbf{Y_o}(R_{\rm match}, E) - \mathbf{Y_i}(R_{\rm match}, E)| = 0$,  (8)

where the subscripts $o$ and $i$ designate, respectively, the logderivative matrices at $R = R_{\rm match}$ obtained from the outward and inward propagations. A good initial estimate of the roots of the determinant in Eq. (8) can be obtained from a prior adiabatic bender calculation.

Equivalently, the energies and wavefunctions can be obtained variationally. Here, each element of the solution matrix, $F_{i,j}(R)$ is expanded in terms of a set of functions

$F_{ij}(R) = \sum_{m=1}^M C_{mi}^{(j)} \chi_m(R) $.  (9)

As in any standard linear variational technique, the energy and eigenfunction of the nth bound state is obtained by diagonalization of the matrix of the full Hamiltonian

$\mathbf{H}(R) = \mathbf{1} \frac{\partial^2}{\partial R^2}+ \mathbf{W}(R)$  (10)

in the $\chi(R)$ basis. The dimensions of this Hamiltonian matrix are $M \times N_{\rm ch}$, where $M$ is the dimensionality of the $\chi(R)$ basis and $N_{\rm ch}$ is the number of internal states (channels). In the present implementation within Hibridon code, the distributed Gaussian basis of Hamilton and Light [6] is used, in which

 $\chi_m(R) = \exp^{-\alpha(R-R_m)^2}$.  (11)

## Coupled States Approximation

Similarly to the treatment of molecular scattering, the [Coupled-States](./Coupled-states.md) approximation can be made for bound states. (This is sometimes called the centrifugal decoupling (CD) approximation [1]). Here, the energies and eigenfunctions refer to a system in which the projection of the total angular momentum along R is assumed to be conserved, so that this projection is a good quantum number.
The Hibridon code will calculate bound-state energies if the flag [`BOUNDC`](./param/boundc.md) is set .TRUE.

Another, powerful program package for the solution of the close-coupled equations for bound states is the [`BOUND`](http://www.dur.ac.uk/~dch0www/Staff/jmh) code developed by J. Hutson.


***

References:
1. J. M. Hutson in Advances in Molecular Vibrations and Collision Dynamics, edited by J. M. Bowman and M. A. Ratner (JAI Press, Greenwich, CT, 1991) 1A, p. 1. 
2. S. L. Holmgren, M. Waldman, and W. Klemperer, J. Chem. Phys. 67, 4414 (1977). 
3. M. H. Alexander, S. Gregurick, and P. J. Dagdigian, J. Chem. Phys. 101, 2887 (1994). 
4. See, for example, J. M. Blatt, J. Comput. Phys. 1, 382 (1967). 
5. D. E. Manolopoulos, Ph. D. thesis, University of Cambridge (1988). 
6. I. P. Hamilton and J. C. Light, J. Chem. Phys. 84, 306 (1986).