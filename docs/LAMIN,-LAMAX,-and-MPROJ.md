In most scattering problems the potential is expanded as

$V(R, \theta) = V_{00}(R) + \sum V_{\lambda}(R) \times f_1(\theta) + \sum V_{\lambda}(R) \times f_2(\theta) +\ ...$


where f_i(\theta) can designate ordinary or associated Legendre polynomials, spherical harmonics, etc.

The number of different sums is given by the parameter NTERM.

For each different sum the minimum value of the expansion index is given by LAMMIN(i), the maximum values of the expansion index, by LAMMAX(i), and the order of the associated Legendre polynomials (or spherical harmonics) by MPROJ(i).

⚠️ The values of NTERM, LAMMIN, and LAMMAX cannot be changed; for any given system they are set in the LOAPOT subroutine.

Likewise, the values of MPROJ cannot be changed; for any given system they are set in the system input subroutine (sy1sg, sy2sg, ...)
