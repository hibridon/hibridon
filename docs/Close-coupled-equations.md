1. [Introduction](#introduction)
2. [Cross sections](#cross-sections)
3. [Wavefunction and Fluxes](#wavefunction-and-fluxes)
4. [Methods](#methods)
5. [References](#references)

## Introduction

The close coupled (CC) equations arise whenever one describes the scattering of systems with internal degrees of freedom.
These equations were introduced into the field of molecular collisions in 1960 by Arthurs and Dalgarno.[1,2]
An excellent description of these equations and the role they play in the quantum theory of inelastic collisions is given in a review article by Secrest.[3] In this Introduction we shall assume that you have a basic familiarity with quantum scattering theory and with the CC equations.

In the close coupled treatment of both scattering and photodissociation, the scattering wavefunction is expanded in a complete set of internal states of the system, usually constructed as direct products of the internal states of one (or both) fragments, multiplied by angular functions which describe the rotation of one collision partner about the other. Let us designate these internal states as (r), where r designates the internal coordinates. Each internal state is called a channel. In the Hibridon code, the channels are labelled by two integer indices, a rotational angular momentum (contained in the array [JOUT](JOUT)) and an additional quantum number or index (contained in the array [INDOUT](INDOUT)).

The full scattering wavefunction is written as:

(1)$$ \Psi(\mathbf{R},\mathbf{r}) = \sum \mathbf{F(R)} \times \Phi(\mathbf{r}) $$

Each column of the $\mathbf{F}(R)$ matrix defines the expansion coefficients for collisions in which the collision partners start out in the particular initial state whose index is that of the selected column.
Substitution of the expansion (1) into the Schrödinger equation, premultiplication by one of the internal states, and integration over $r$ gives rise to a set of coupled ordinary differential equations for the expansion coefficients $\mathbf{F}(R)$.

The general structure of these coupled second-order differential equations is expressed by the matrix equation:

(2)$$ \left[ \mathbf{1}\frac{\partial^2}{\partial R²} + \mathbf{W}(R)\right] \mathbf{F}(R) = 0 $$

Here $\mathbf{1}$ designates the identity matrix, $R$ is the interparticle distance, and the matrix $\mathbf{W}(R)$ is given by:

(3)$$ \mathbf{W}(R)= \mathbf{k}^2 - \mathbf{l}^2 - \left( \frac{2m}{\hbar^2} \right)  \mathbf{V}(R)$$


where $\hbar$ designates Planck's constant divided by $2\pi$, $m$ is the reduced mass of the collision system and $\mathbf{k}^2$ and $\mathbf{l}^2$ designate, respectively, the (diagonal) matrices of the wavevector and the relative orbital angular momentum of the collision partners. We have

(4)$$ (k^2)_{ii} = \left( \frac{2m}{\hbar^2} \right) (E-e_i) $$

where $E$ is the total energy and $e_i$ is the internal energy of the $i$<sup>th</sup> channel. Also, we have

(5)$$ (l^2)_{ii} = \left( \frac{\hbar^2}{2mR^2} \right) l_i(l_i + 1) $$

where $l_i$ is the relative orbital angular momentum in the $i$<sup>th</sup> channel. In Eq. (3) the matrix $\mathbf{V}(R)$ is the (full, symmetric) matrix of the coupling potential.

Diagonalization of the $\mathbf{W}(R)$ matrix yields the diagonal matrix of adiabatic wavevectors $\mathbf{k}(R)$. The eigenvectors define the locally adiabatic states, which are transformations of the internal states used to expand the scattering wavefunction, (r). If $\mathbf{C}(R)$ designates the matrix of eigenvectors, column ordered, then the diagonal matrix of adiabatic energies is defined as

(6)$$ \mathbf{e}(R) = \mathbf{C}(R)^T\mathbf{V}(R)\mathbf{C}(R) $$

One obtains numerically the matrix of solutions $\mathbf{F}(R)$ by outward propagation. You start this propagation at a value of the interparticle distance $R$ = $R_{start}$ which lies well inside the innermost classical turning point. Once you have propagated $\mathbf{F}(R)$ out to a value of $R$ which is so large that the potential $\mathbf{V}(R)$ is negligible, compared to the wavevectors $\mathbf{k}^2$, you can then match $\mathbf{F}(R)$ to the known [asymptotic form](Boundary) and obtain the S matrix. This has to be done over and over at many values of the total angular momentum $J_{tot}$. In a semiclassical description the total angular momentum corresponds to the impact parameter b. From the $\mathbf{S}$ matrix at all these values of $J_{tot}$, you can calculate differential and integral cross sections.[2,3]

## Cross sections

Within the CC formulation, the integral cross section for a transition from an initial state $i$ to a final state $f$ is given by

(7)$$\sigma_{i \rightarrow f} = \frac{\pi}{(2j_i + 1)k_i^2} \sum_{J_{tot}, l, p} \left|T_{ij}^{J_{tot}} \right|^2$$
 
where $k_i$ is the wavevector of the initial state and the sum runs over:

- all values of the total angular momentum for which the $\mathbf{S}$ matrix elements differ from unity. You must ensure that the parameter [JTOT2](JTOT) has been set large enough so that increasing its value will not significantly change the cross section(s) of interest,
- all values of the orbital angular momentum $l$ allowed by the triangular rule
- both values ($p$ = +1 and -1) of the total parity of the scattering wavefunction. The total parity is related to the important input parameter [JLPAR](JLPAR). Note that for full close-coupled determinations of either integral or differential cross sections, the calculations must be carried out for both values of [JLPAR](JLPAR).

Here, the $\mathbf{T}$ or transition matrix is defined as:

(8)$$ \mathbf{T} = \mathbf{1} - \mathbf{S} $$

where $\mathbf{1}$ is the unit matrix. At large $J_{tot}$ the $\mathbf{T}$ matrix goes to zero as the centrifugal potential becomes so large that the colliding particles are kept beyond the range of the interaction potential. This defines the range of total angular momentum for which scattering calculations need be done. The minimum and maximum values of the total angular momentum for which the calculation is done are set by the parameters [JTOT1](JTOT) and [JTOT2](JTOT), respectively.

See Refs. 1 and 2 for an expression for the differential cross section equivalent to Eq. (7).

In general, the CC equations are block diagonal in the overall parity of the scattering wavefunction. To obtain integral and/or differential cross sections it is necessary to carry out calculations for both values of this parity (this is ensured by setting [JLPAR](JLPAR)=0; see the [JLPAR](JLPAR) page for more information).

Equation (7) can be written equivalently in terms of partial cross sections

(9)$$\sigma_{i \rightarrow f} = \sum_{J_{tot}} \sigma_{i \rightarrow f}^{J_{tot}}$$

where the partial cross sections, which can be calculated with the command [PARTC](PARTC), are defined by:

(10)$$\sigma_{i \rightarrow f}^{J_{tot}} = \frac{\pi}{(2j_i + 1)k_i^2} \sum_{l, p} \left|T_{ij}^{J_{tot}} \right|^2$$

In a semiclassical formulation, the integral cross section is written as

(11)$$\sigma_{i \rightarrow f} = \int_0^\infty 2\pi b P_{i \rightarrow f}(b)db$$


in terms of a transition probability which depends on the impact parameter $b$. The partial cross section is thus equivalent to this semi-classical transition probability, as follows:

(12)$$\sigma_{i \rightarrow f}^{J_{tot}} = 2\pi b P_{i \rightarrow f}(b)db$$

## Wavefunction and Fluxes

In addition to the determination of cross sections, which depend on the asymptotic form of the expansion coefficients $\mathbf{F}(R)$ in Eq. (1), the Hibridon code allows you determine the R dependence of these coefficients. In addition, one can also determine the flux associated with the scattering and photodissociation wavefunctions, which is a more meaningful physical quantity. For a further discussion see Ref. 4 and 5.)

## Methods

The determination of differential and/or integral cross sections involves three steps:
The development of subroutines to calculate the potential matrix $\mathbf{V}(R)$ for a particular collision system.
Solution of the CC equations to obtain and store the $\mathbf{S}$ matrix elements, at both values of the parity index [JLPAR](JLPAR).
The subsequent calculation of differential and/or integral cross sections for the transitions which interest you.
The major part of this manual is devoted to the description of a complex family of subroutines for solution of the CC equations - step 2 above.
   Historically, there have been many algorithms developed to solve these equations. These algorithms can be grouped into two categories[3]:

- Solution-following methods. In these methods you approximate the matrix of solutions $\mathbf{F}(R)$ by a power series and then solve Eq. (1) exactly. This is similar in spirit to the usual numerical techniques for solution of ordinary differential equations (Runge-Kutta, Euler, Adams-Moulton).
- Potential-following methods. In these methods the matrix $\mathbf{V}(R)$ is approximated by a sequence of constant or linear segments. In these local regions the approximated CC equations can be solved exactly.

In solution-following methods the solution is approximated while the potential is retained exactly. On the other hand, in potential-following methods the potential is approximated but the solution (to this approximate potential) is exact.

   No one method is superior at all values of the internuclear separation. Rather, it is best to combine a solution-following method at short-range ($R$ small), where most intermolecular potentials vary rapidly, with a potential-following method at longer range, where the potential varies more slowly but where for many problems the solution can be highly oscillatory. This combination of two methods is called a hybrid integrator.

   The Hibridon program package uses a particular hybrid integrator. The solution following method used at short range is based on the log-derivative propagator of Johnson,[6-8] as modified recently by Manolopoulos.[9] This propagator is designated LOGD. The potential-following method used at long-range is based on the linear-reference potential of Gordon,[10,11] as modified recently by Alexander and Manolopoulos.[12,13] This propagator is designated AIRY. The Hibridon code combines these two fast algorithms (LOGD and AIRY). Both are fast and exceptionally stable. To a large degree the numerical stability is obtained by propagation of the logarithmic derivative of the solution matrix $\mathbf{F}(R)$, namely

(13)$$ \mathbf{Y}(R) = \mathbf{F}'(R) \mathbf{F}(R)^{-1}$$

rather than the solution matrix itself.

   Another, powerful program package for the solution of the close coupled equations is the [MOLSCAT](http://www.giss.nasa.gov/tools/molscat/) code developed by S. Green and maintained by J. Hutson.

## References
1. A. Arthurs and A. Dalgarno, Proc. Roy. Soc. (London Ser.) **A256**, 540 (1960). 
2. W. A. Lester, Jr., Meth. Comput. Phys. **10**, 211 (1971). 
3. D. Secrest, in Atom-Molecule Collision Theory: A Guide for the Experimentalist, edited by R. B. Bernstein (Plenum, New York, l979) p. 265. 
4. M. H. Alexander, J. Chem. Phys. **95**, 8931 (1991); **96**, 6672 (1992). 
5. D. E. Manolopoulos and M. H. Alexander, J. Chem. Phys. **97**, 2527 (1992); M. H. Alexander, C. Rist, and D. E. Manolopoulos, J. Chem. Phys. **97**, 4836 (1992). 
6. B. R. Johnson, J. Comput. Phys. **13**, 445 (1973). 
7. B. R. Johnson, Proceedings of the NRCC Workshop on Algorithms and Computer Codes in Atomic and Molecular Scattering Theory, edited by L. D. Thomas (Lawrence Berkeley Laboratory, CA Report LBL-9501 l979) pp. 86-92 (Vol. I) and p. 52 (Vol. II). 
8. F. Mrugala and D. Secrest, J. Chem. Phys. **78**, 5954 (1983). 
9. D. E. Manolopoulos, J. Chem. Phys. **85**, 6425 (1986). 
10. R. G. Gordon, Meth. Comput. Phys. **10**, 81 (1971). 
11. R. G. Gordon, J. Chem. Phys. 51, **14** (1969). 
12. M. H. Alexander, J. Chem. Phys. **81**, 4510 (1984). 
13. M. H. Alexander and D. E. Manolopoulos, J. Chem. Phys. **86**, 2044 (1987).