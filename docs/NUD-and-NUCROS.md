NUD is the step size for the angular momentum projection in a [Coupled States](Coupled-States) calculation. Scattering calculations are performed at every NUD value of the projection quantum number nu lying between (and including) [NUMIN](NUMIN-and-NUMAX) and [NUMAX](NUMIN-and-NUMAX). For collisions involving a molecule with half-integer spin ([flag](Flags) [FLAGHF](FLAGHF)=.TRUE.) then the total angular momentum and its projection is half-integer (although the parameter NU remains integer but is now equal to the projection minus 1/2). In this case the values of the projection quantum number range from NUMIN+1/2 to NUMAX+1/2 in steps of NUD.

If the [flags](Flags) [WRXSEC](PRXSEC-and-WRXSEC) and/or [PRXSEC](PRXSEC-and-WRXSEC) are .TRUE., then the program accumulates the contribution to the integral cross sections due to all partial waves lying between [JTOT1](JTOT) and [JTOT2](JTOT). If [JTOTD](JTOT)=1, this is done by adding the degeneracy weighted sums of the squares of the T matrix at each partial wave.

In a [CS](Coupled-States) calculation the integral cross sections can be determined in two ways:

If NUCROS = .FALSE., at each partial wave (for each value of $\bar{L}$ ), the squares of the T matrix for each value of the projection index are summed to give a partial cross section indexed only in $\bar{L}$. 

Then, similarly to a [close-coupled](Close-coupled-equations) calculation the program accumulates the contribution to the integral cross sections due to all partial waves lying between [JTOT1](JTOT) and [JTOT2](JTOT).

If NUCROS = .TRUE., the summation over partial wave is done before the summation over the values of the projection index. The first summation yields partial cross sections indexed only in the projection index. These are then summed to give the integral cross section.

In both cases if

  - NUD = 1, the sum over projection quantum number $\bar{L}$ is done by adding the contribution from each partial wave.
  - NUD > 1, the sum over partial wave $\bar{L}$ is done using a Simpson's rule interpolation of these NU-indexed integral cross sections.
